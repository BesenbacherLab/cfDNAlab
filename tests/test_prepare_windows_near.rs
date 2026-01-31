//! TODO: These tests are completely unvalidated. Go through each and every one manually.

#![cfg(feature = "cmd_prepare_windows")]

mod tests_prepare_windows_near {
    use anyhow::Result;
    use cfdnalab::commands::prepare_windows::chunk::apply_near_annotations;
    use cfdnalab::commands::prepare_windows::config::{
        CoordinateSet, DistSign, NearDirection, NearEdge, NearTiePolicy, PrepareConfig,
    };
    use cfdnalab::commands::prepare_windows::labels::{
        LabelSchema, LabelTuple, build_tuple_compositions, render_label_for_key,
    };
    use cfdnalab::commands::prepare_windows::near_file::NearDuplicatesPolicy;
    use cfdnalab::commands::prepare_windows::near_file::{
        NearChrom, NearInterval, NearWindowSide, NearestDistance, NearestResult, Strand,
        nearest_edge_distance,
    };
    use cfdnalab::commands::prepare_windows::parsers::parse_distance_bins;
    use cfdnalab::commands::prepare_windows::prepare_windows::Window;
    use std::io::Write;
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn build_window(chrom: &str, start: u32, end: u32, input_label: &str) -> Window {
        Window {
            chrom: Arc::from(chrom),
            original_start: start,
            original_end: end,
            resized_start: start,
            resized_end: end,
            merged: false,
            label_tuples: vec![LabelTuple::new(input_label.to_string())],
            group_key: input_label.to_string(),
            score: None,
        }
    }

    fn label_schema_from_compose(tokens: &[&str]) -> LabelSchema {
        let specs = tokens
            .iter()
            .map(|spec| spec.parse().expect("compose spec"))
            .collect::<Vec<_>>();
        LabelSchema::new(&specs).expect("schema")
    }

    fn make_near_index(
        intervals: Vec<NearInterval>,
    ) -> cfdnalab::commands::prepare_windows::near_file::NearIndex {
        let mut idx = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        idx.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals,
                cursor: 0,
            },
        );
        idx
    }

    #[test]
    fn nearest_edge_distance_reports_negative_upstream_on_plus_strand() {
        // Arrange
        // Interval on + strand: upstream edge is left coordinate 100
        // Window lies to the left (50-60), closest point to edge 100 is end=60 -> distance -40
        let interval = NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            50,
            60,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, -40);
            assert_eq!(window_side, NearWindowSide::Upstream);
        } else {
            panic!("expected single upstream hit");
        }
    }

    #[test]
    fn nearest_edge_distance_reports_negative_upstream_on_minus_strand() {
        // Arrange
        // Interval on - strand: upstream edge is right coordinate 110
        // Window lies to the right (130-140); genomic distance is +20, but strand flips sign
        let interval = NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Minus,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            130,
            140,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, -20); // flipped by strand
            assert_eq!(window_side, NearWindowSide::Upstream); // genomic downstream becomes upstream on -
        } else {
            panic!("expected single upstream hit");
        }
    }

    #[test]
    fn nearest_edge_distance_returns_tie_when_equidistant() {
        // Arrange
        // Window 15-25 is 5 bp away from right edge of left interval (10) and left edge of right interval (30)
        let intervals = vec![
            NearInterval {
                start: 0,
                end: 10,
                group_id: None,
                strand: Strand::Plus,
            },
            NearInterval {
                start: 30,
                end: 40,
                group_id: None,
                strand: Strand::Plus,
            },
        ];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            15,
            25,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert
        match result {
            NearestResult::Tie(tie) => {
                let upstream = tie.left.expect("left");
                let downstream = tie.right.expect("right");
                assert_eq!(upstream.distance, -5);
                assert_eq!(upstream.window_side, NearWindowSide::Upstream); // window lies upstream of right interval
                assert_eq!(downstream.distance, 5);
                assert_eq!(downstream.window_side, NearWindowSide::Downstream); // window lies downstream of left interval
            }
            _ => panic!("expected tie"),
        }
    }

    #[test]
    fn apply_near_annotations_sets_no_near_labels_when_bins_and_no_distance_max() -> Result<()> {
        // Arrange
        // No near intervals on any chromosome, but bins are requested
        // Should retain window and label with [NONE]/[NO-NEAR]
        let windows = vec![build_window("chr1", 10, 20, "input1")];
        let mut near_index = Some(Default::default());
        let mut cfg = PrepareConfig::default();
        cfg.near_group_cols = vec!["3".to_string()]; // triggers near-name labeling path
        let bins = parse_distance_bins(&vec!["prox:<100".to_string()])?;

        // Act
        let result = apply_near_annotations(
            windows,
            &mut near_index,
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        assert_eq!(tuple.near_side.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.near_name.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.bin.as_deref(), Some("[NO-NEAR]"));
        Ok(())
    }

    #[test]
    fn apply_near_annotations_drops_when_direction_mismatch_and_distance_max() {
        // Arrange
        // Near interval exists downstream, but we only accept upstream; distance_max forces drop
        let interval = NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        };
        let mut near_index = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        near_index.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals: vec![interval],
                cursor: 0,
            },
        );

        let windows = vec![build_window("chr1", 120, 130, "input1")];
        let mut cfg = PrepareConfig::default();
        cfg.distance_max = Some(50);
        cfg.near_direction = NearDirection::Upstream; // window is downstream, so no match

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert!(result.is_empty());
    }

    #[test]
    fn apply_near_annotations_builds_composed_label_for_tie_annotations() -> Result<()> {
        // Arrange
        // Window equidistant to upstream and downstream sites with group names; tie annotate
        // produces two tuples that render to comma-joined composition
        let intervals = vec![
            NearInterval {
                start: 0,
                end: 10,
                group_id: Some(0),
                strand: Strand::Plus,
            },
            NearInterval {
                start: 30,
                end: 40,
                group_id: Some(1),
                strand: Strand::Plus,
            },
        ];
        let mut near_index = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        near_index.group_id_to_name = vec!["UP".to_string(), "DN".to_string()];
        near_index.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals,
                cursor: 0,
            },
        );

        let windows = vec![build_window("chr1", 15, 25, "input1")];
        let mut cfg = PrepareConfig::default();
        cfg.near_ties = NearTiePolicy::Annotate;
        cfg.distance_sign = DistSign::Signed;
        cfg.out_labels = vec!["near".to_string()];
        cfg.compose = vec!["near=win-direction,near-name".parse().expect("compose")];

        let schema = label_schema_from_compose(&["near=win-direction,near-name"]);
        let near_key = schema.resolve_key("near")?;
        let out_keys = schema.resolve_keys(&cfg.out_labels)?;

        // Act
        let mut annotated = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Render composed label manually using same helper as pipeline
        for window in annotated.iter_mut() {
            let compositions = build_tuple_compositions(&window.label_tuples, &schema);
            window.group_key =
                render_label_for_key(&window.label_tuples, &compositions, &near_key, &schema);
            // project to output labels
            let rendered =
                render_label_for_key(&window.label_tuples, &compositions, &out_keys[0], &schema);
            window.group_key = rendered;
        }

        // Assert
        // Upstream interval appears first; downstream interval second
        assert_eq!(annotated.len(), 1);
        assert_eq!(annotated[0].group_key, "-.DN,+.UP");
        Ok(())
    }

    #[test]
    fn nearest_edge_distance_reports_overlap_as_zero() {
        // Arrange
        // Window overlaps interval (5-15 vs 10-20) so distance must be 0 regardless of edge/direction
        let interval = NearInterval {
            start: 10,
            end: 20,
            group_id: None,
            strand: Strand::Plus,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            5,
            15,
            &mut chrom,
            &NearEdge::Right,
            &NearDirection::Upstream,
            false,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, 0);
            assert_eq!(window_side, NearWindowSide::Overlap);
        } else {
            panic!("expected overlap hit");
        }
    }

    #[test]
    fn nearest_edge_distance_falls_back_for_unknown_strand_on_directional_edge() {
        // Arrange
        // Unknown strand causes upstream edge mode to behave like nearest; window right of interval so side downstream
        let interval = NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Unknown,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            130,
            140,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            // Falls back to nearest edge: right edge at 110 is closest, window downstream so +20
            assert_eq!(distance, 20);
            assert_eq!(window_side, NearWindowSide::Downstream);
        } else {
            panic!("expected downstream hit");
        }
    }

    #[test]
    fn apply_near_annotations_uses_resized_coordinates_for_distance_binning() -> Result<()> {
        // Arrange
        // Resized window 15-35 is 5 bp from interval edge (10), original 20-30 would be 10 bp
        // distance_from=resized should pick bin based on 5 bp
        let mut near_index = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        near_index.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals: vec![NearInterval {
                    start: 10,
                    end: 12,
                    group_id: None,
                    strand: Strand::Plus,
                }],
                cursor: 0,
            },
        );

        let mut cfg = PrepareConfig::default();
        cfg.distance_from = CoordinateSet::Resized;
        cfg.distance_bins = Some(vec!["near:<7".to_string(), "far:>=7".to_string()]);
        cfg.flank = Some(vec![5, 5]); // resize via flank for test fixture consistency
        cfg.oob = cfdnalab::commands::prepare_windows::config::OobPolicy::Allow;
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        let windows = vec![build_window("chr1", 20, 30, "A")]; // resized to 15-35

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        // 5 bp distance should fall in "near"
        assert_eq!(tuple.bin.as_deref(), Some("near"));
        Ok(())
    }

    #[test]
    fn apply_near_annotations_uses_original_coordinates_for_distance_binning() -> Result<()> {
        // Arrange
        // Original window 20-30 is 10 bp from interval edge (10), resized 15-35 would be 5 bp
        // distance_from=original should pick bin based on 10 bp
        let mut near_index = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        near_index.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals: vec![NearInterval {
                    start: 10,
                    end: 12,
                    group_id: None,
                    strand: Strand::Plus,
                }],
                cursor: 0,
            },
        );

        let mut cfg = PrepareConfig::default();
        cfg.distance_from = CoordinateSet::Original;
        cfg.distance_bins = Some(vec!["near:<7".to_string(), "far:>=7".to_string()]);
        cfg.flank = Some(vec![5, 5]);
        cfg.oob = cfdnalab::commands::prepare_windows::config::OobPolicy::Allow;
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        let windows = vec![build_window("chr1", 20, 30, "A")]; // resized to 15-35 but binning uses 20-30

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        // 10 bp distance should fall in "far"
        assert_eq!(tuple.bin.as_deref(), Some("far"));
        Ok(())
    }

    #[test]
    fn apply_near_annotations_fills_missing_near_group_with_na() -> Result<()> {
        // Arrange
        // Near interval with empty group column should emit near-name [NA]
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t0\t10\t")?;

        let index = cfdnalab::commands::prepare_windows::near_file::load_near_index(
            file.path(),
            '\t',
            false,
            None,
            Some(&[3]),
            false,
            cfdnalab::commands::prepare_windows::near_file::NearDuplicatesPolicy::Error,
        )?;

        let windows = vec![build_window("chr1", 20, 22, "A")];
        let mut cfg = PrepareConfig::default();
        cfg.near_group_cols = vec!["3".to_string()];

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        assert_eq!(tuple.near_name.as_deref(), Some("[NA]"));
        Ok(())
    }

    #[test]
    fn load_near_index_merges_duplicates_with_strand_considered() -> Result<()> {
        // Arrange
        // Two identical intervals with same strand should merge groups when policy merge and strand considered
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t10\t20\t+\tA")?;
        writeln!(file, "chr1\t10\t20\t+\tB")?;

        // Act
        let index = cfdnalab::commands::prepare_windows::near_file::load_near_index(
            file.path(),
            '\t',
            false,
            Some(3),
            Some(&[4]),
            true, // consider strand for upstream/downstream edge modes
            cfdnalab::commands::prepare_windows::near_file::NearDuplicatesPolicy::Merge,
        )?;

        // Assert
        let chr1 = index.per_chrom.get("chr1").expect("chr1");
        assert_eq!(chr1.intervals.len(), 1);
        let gid = chr1.intervals[0].group_id.expect("group id");
        let name = &index.group_id_to_name[gid as usize];
        assert_eq!(name, "A__B");
        Ok(())
    }

    #[test]
    fn load_near_index_keep_first_drops_subsequent_duplicates() -> Result<()> {
        // Arrange: two identical intervals; keep-first retains the first only
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t10\t20\t+\tX")?;
        writeln!(file, "chr1\t10\t20\t+\tY")?;

        // Act
        let index = cfdnalab::commands::prepare_windows::near_file::load_near_index(
            file.path(),
            '\t',
            false,
            Some(3),
            Some(&[4]),
            true,
            NearDuplicatesPolicy::KeepFirst,
        )?;

        // Assert: only one interval remains and group is from first record
        let chr1 = index.per_chrom.get("chr1").expect("chr1");
        assert_eq!(chr1.intervals.len(), 1);
        let gid = chr1.intervals[0].group_id.expect("group id");
        let name = &index.group_id_to_name[gid as usize];
        assert_eq!(name, "X");
        Ok(())
    }

    #[test]
    fn load_near_index_drop_all_removes_duplicate_run() -> Result<()> {
        // Arrange: two identical intervals; drop-all removes the entire run
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t10\t20\t+\tX")?;
        writeln!(file, "chr1\t10\t20\t+\tY")?;

        // Act
        let index = cfdnalab::commands::prepare_windows::near_file::load_near_index(
            file.path(),
            '\t',
            false,
            Some(3),
            Some(&[4]),
            true,
            NearDuplicatesPolicy::DropAll,
        )?;

        // Assert: no intervals remain for chr1
        let chr1 = index.per_chrom.get("chr1").expect("chr1");
        assert!(chr1.intervals.is_empty());
        Ok(())
    }

    #[test]
    fn near_direction_upstream_filters_out_downstream_hit() {
        // Arrange: window downstream of near interval, direction=upstream should block it
        let intervals = vec![NearInterval {
            start: 0,
            end: 10,
            group_id: None,
            strand: Strand::Plus,
        }];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            20,
            25,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Upstream,
            true,
        );

        // Assert
        assert!(result.is_none());
    }

    #[test]
    fn near_direction_downstream_filters_out_upstream_hit() {
        // Arrange: window upstream of near interval, direction=downstream should block it
        let intervals = vec![NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        }];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            50,
            60,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Downstream,
            true,
        );

        // Assert
        assert!(result.is_none());
    }

    #[test]
    fn bins_applied_when_no_hit_sets_no_near_labels() -> Result<()> {
        // Arrange: chromosome present in near index but window falls outside direction/edge filter, distance_max unset
        let intervals = vec![NearInterval {
            start: 0,
            end: 10,
            group_id: None,
            strand: Strand::Plus,
        }];
        let near_index = make_near_index(intervals);
        let mut cfg = PrepareConfig::default();
        cfg.near_direction = NearDirection::Upstream; // window downstream
        cfg.distance_bins = Some(vec!["prox:<5".to_string()]);
        cfg.near_group_cols = vec!["3".to_string()]; // force near-name path

        let windows = vec![build_window("chr1", 50, 60, "A")];
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert: kept window with [NONE]/[NO-NEAR] because distance_max is None
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        assert_eq!(tuple.near_side.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.near_name.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.bin.as_deref(), Some("[NO-NEAR]"));
        Ok(())
    }

    #[test]
    fn touch_is_treated_as_overlap_with_zero_distance() {
        // Arrange: window ends exactly where interval starts; treated as overlap with distance 0
        let intervals = vec![NearInterval {
            start: 20,
            end: 30,
            group_id: None,
            strand: Strand::Plus,
        }];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            10,
            20,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, 0);
            assert_eq!(window_side, NearWindowSide::Overlap);
        } else {
            panic!("expected touch -> overlap hit");
        }
    }

    #[test]
    fn upstream_edge_respects_strand_when_strand_column_present() {
        // Arrange: same genomic interval annotated on + and -; upstream edge differs by strand
        let intervals = vec![
            NearInterval {
                start: 100,
                end: 110,
                group_id: None,
                strand: Strand::Plus,
            },
            NearInterval {
                start: 200,
                end: 210,
                group_id: None,
                strand: Strand::Minus,
            },
        ];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act: window between them (120-125) should be downstream of first (+ strand), upstream of second (- strand flips)
        let first = nearest_edge_distance(
            90,
            95,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit1");
        // Move cursor forward and query near the second interval
        let second = nearest_edge_distance(
            220,
            225,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit2");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = first
        {
            assert!(distance > 0); // window 90-95 is downstream of upstream edge 100 on +
            assert_eq!(window_side, NearWindowSide::Downstream);
        } else {
            panic!("expected downstream hit vs + strand upstream edge");
        }
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = second
        {
            // On - strand upstream edge is right coordinate 210; window 220-225 is upstream (positive genomic) but flips sign
            assert!(distance < 0);
            assert_eq!(window_side, NearWindowSide::Upstream);
        } else {
            panic!("expected upstream hit vs - strand upstream edge");
        }
    }

    #[test]
    fn downstream_edge_respects_strand_when_strand_column_present() {
        // Arrange: same as prior but downstream edge selection
        let intervals = vec![
            NearInterval {
                start: 100,
                end: 110,
                group_id: None,
                strand: Strand::Plus,
            },
            NearInterval {
                start: 200,
                end: 210,
                group_id: None,
                strand: Strand::Minus,
            },
        ];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act: window left of first should be upstream; window right of second should be downstream after strand flip
        let first = nearest_edge_distance(
            90,
            95,
            &mut chrom,
            &NearEdge::Downstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit1");
        let second = nearest_edge_distance(
            220,
            225,
            &mut chrom,
            &NearEdge::Downstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit2");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = first
        {
            assert!(distance < 0); // downstream edge on + is right edge 110; window left -> upstream negative
            assert_eq!(window_side, NearWindowSide::Upstream);
        } else {
            panic!("expected upstream hit vs + strand downstream edge");
        }
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = second
        {
            // On - strand downstream edge is left coordinate 200; window right -> downstream positive
            assert!(distance > 0);
            assert_eq!(window_side, NearWindowSide::Downstream);
        } else {
            panic!("expected downstream hit vs - strand downstream edge");
        }
    }

    #[test]
    fn distance_max_with_overlap_keeps_overlap() {
        // Arrange: overlap distance is 0, so distance_max should keep it even when set small
        let intervals = vec![NearInterval {
            start: 10,
            end: 20,
            group_id: None,
            strand: Strand::Plus,
        }];
        let near_index = make_near_index(intervals);
        let mut cfg = PrepareConfig::default();
        cfg.distance_max = Some(1);

        let windows = vec![build_window("chr1", 15, 18, "A")]; // overlap

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn unknown_strand_upstream_falls_back_to_nearest_for_distance() {
        // Arrange: unknown strand and upstream edge should behave like nearest
        let intervals = vec![NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Unknown,
        }];
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            80,
            82,
            &mut chrom,
            &NearEdge::Upstream,
            &NearDirection::Both,
            true,
        )
        .expect("hit");

        // Assert: nearest edge is start=100, window is upstream so negative distance
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert!(distance < 0);
            assert_eq!(window_side, NearWindowSide::Upstream);
        } else {
            panic!("expected upstream hit with unknown strand");
        }
    }

    #[test]
    fn large_interval_list_cursor_stays_linear() {
        // Arrange: many intervals; ensure cursor advances without quadratic behavior by checking it reaches near end
        let mut intervals: Vec<NearInterval> = Vec::new();
        for i in 0..1000 {
            intervals.push(NearInterval {
                start: i * 100,
                end: i * 100 + 10,
                group_id: None,
                strand: Strand::Plus,
            });
        }
        let mut chrom = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act: query far downstream so cursor should move to last interval
        let _ = nearest_edge_distance(
            100_000,
            100_010,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            false,
        );

        // Assert: cursor at or near last element (999)
        assert!(chrom.cursor >= 998);
    }

    #[test]
    fn label_rendering_composes_single_value_when_only_input_differs() -> Result<()> {
        // Arrange
        // Two tuples differ only in input; composition depends on win-direction so it should collapse to one value
        let tuples = vec![
            LabelTuple {
                input: "A".into(),
                near_side: Some("-".into()),
                near_name: Some("G".into()),
                bin: None,
                cluster: None,
            },
            LabelTuple {
                input: "B".into(),
                near_side: Some("-".into()),
                near_name: Some("G".into()),
                bin: None,
                cluster: None,
            },
        ];
        let schema = label_schema_from_compose(&["core=win-direction,near-name"]);
        let core_key = schema.resolve_key("core")?;
        let compositions = build_tuple_compositions(&tuples, &schema);

        // Act
        let rendered = render_label_for_key(&tuples, &compositions, &core_key, &schema);

        // Assert
        // Only input differs; composition does not depend on input, so emit single value "-.G"
        assert_eq!(rendered, "-.G");
        Ok(())
    }

    #[test]
    fn label_rendering_keeps_order_when_parts_differ() -> Result<()> {
        // Arrange
        // Tuples differ in win-direction; composition uses win-direction+near-name so it should comma-join in tuple order
        let tuples = vec![
            LabelTuple {
                input: "A".into(),
                near_side: Some("-".into()),
                near_name: Some("UP".into()),
                bin: None,
                cluster: None,
            },
            LabelTuple {
                input: "A".into(),
                near_side: Some("+".into()),
                near_name: Some("DN".into()),
                bin: None,
                cluster: None,
            },
        ];
        let schema = label_schema_from_compose(&["near=win-direction,near-name"]);
        let near_key = schema.resolve_key("near")?;
        let compositions = build_tuple_compositions(&tuples, &schema);

        // Act
        let rendered = render_label_for_key(&tuples, &compositions, &near_key, &schema);

        // Assert
        assert_eq!(rendered, "-.UP,+.DN");
        Ok(())
    }

    #[test]
    fn cursor_rewinds_when_window_moves_backward() {
        // Arrange
        // First call moves cursor forward; second call with earlier window should still find upstream interval
        let interval = NearInterval {
            start: 10,
            end: 20,
            group_id: None,
            strand: Strand::Plus,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Prime cursor with forward window
        let _ = nearest_edge_distance(
            30,
            40,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            false,
        );

        // Act: earlier window should still hit
        let result = nearest_edge_distance(
            5,
            8,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Both,
            false,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, 2); // nearest edge at 10 vs window end 8 -> distance 2
            assert_eq!(window_side, NearWindowSide::Upstream);
        } else {
            panic!("expected upstream hit after cursor rewind");
        }
    }

    #[test]
    fn overlapping_window_kept_when_direction_is_upstream_only() {
        // Arrange
        // Overlap should always count as a hit even if direction filter would otherwise block
        let interval = NearInterval {
            start: 10,
            end: 20,
            group_id: None,
            strand: Strand::Plus,
        };
        let mut chrom = NearChrom {
            intervals: vec![interval],
            cursor: 0,
        };

        // Act
        let result = nearest_edge_distance(
            15,
            25,
            &mut chrom,
            &NearEdge::Nearest,
            &NearDirection::Upstream,
            true,
        )
        .expect("hit");

        // Assert
        if let NearestResult::Single(NearestDistance {
            distance,
            window_side,
            ..
        }) = result
        {
            assert_eq!(distance, 0);
            assert_eq!(window_side, NearWindowSide::Overlap);
        } else {
            panic!("expected overlap hit");
        }
    }

    #[test]
    fn distance_max_drops_when_direction_and_edge_block_hits() {
        // Arrange: interval exists but edge/direction filter makes it invisible; distance_max should drop window
        let intervals = vec![NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        }];
        let near_index = make_near_index(intervals);
        let mut cfg = PrepareConfig::default();
        cfg.distance_max = Some(50);
        cfg.near_direction = NearDirection::Upstream;
        cfg.near_edge = NearEdge::Downstream; // window upstream, edge choice blocks

        let windows = vec![build_window("chr1", 10, 20, "A")]; // upstream of interval

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert!(result.is_empty());
    }

    #[test]
    fn cross_chromosome_no_near_sets_none_labels() -> Result<()> {
        // Arrange: chr1 has intervals; chr2 lacks them. Chr2 windows should be labeled [NONE]/[NO-NEAR] with bins
        let mut idx = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        idx.per_chrom.insert(
            "chr1".to_string(),
            NearChrom {
                intervals: vec![NearInterval {
                    start: 0,
                    end: 10,
                    group_id: None,
                    strand: Strand::Plus,
                }],
                cursor: 0,
            },
        );

        let mut cfg = PrepareConfig::default();
        cfg.near_group_cols = vec!["3".to_string()];
        cfg.distance_bins = Some(vec!["prox:<5".to_string()]);
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        let windows = vec![build_window("chr2", 50, 60, "A")];

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(idx),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        assert_eq!(tuple.near_side.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.near_name.as_deref(), Some("[NONE]"));
        assert_eq!(tuple.bin.as_deref(), Some("[NO-NEAR]"));
        Ok(())
    }

    #[test]
    fn chromosome_with_no_near_and_distance_max_drops_window() {
        // Arrange: chromosome missing from near index; distance_max should drop window
        let idx = cfdnalab::commands::prepare_windows::near_file::NearIndex::default();
        let mut cfg = PrepareConfig::default();
        cfg.distance_max = Some(100);
        let windows = vec![build_window("chrMissing", 0, 10, "A")];

        // Act
        let result =
            apply_near_annotations(windows, &mut Some(idx), &cfg, None, CoordinateSet::Resized);

        // Assert
        assert!(result.is_empty());
    }

    #[test]
    fn signed_bins_split_around_zero() -> Result<()> {
        // Arrange: signed mode, bins straddle zero to classify upstream/overlap/downstream
        let near_index = make_near_index(vec![NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        }]);
        let mut cfg = PrepareConfig::default();
        cfg.distance_sign = DistSign::Signed;
        cfg.distance_bins = Some(vec![
            "up:<0".to_string(),
            "at:0-0".to_string(),
            "down:>0".to_string(),
        ]);
        cfg.out_labels = vec!["win-direction".to_string(), "bin".to_string()];
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        // Upstream window
        let upstream = vec![build_window("chr1", 50, 60, "A")];
        let upstream_res = apply_near_annotations(
            upstream,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(upstream_res[0].label_tuples[0].bin.as_deref(), Some("up"));

        // Overlap window
        let overlap = vec![build_window("chr1", 105, 108, "A")];
        let overlap_res = apply_near_annotations(
            overlap,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(overlap_res[0].label_tuples[0].bin.as_deref(), Some("at"));

        // Downstream window
        let downstream = vec![build_window("chr1", 150, 160, "A")];
        let downstream_res = apply_near_annotations(
            downstream,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            downstream_res[0].label_tuples[0].bin.as_deref(),
            Some("down")
        );
        Ok(())
    }

    #[test]
    fn absolute_mode_still_emits_direction_prefix() -> Result<()> {
        // Arrange: distance_sign absolute but win-direction should include +/- based on position
        let near_index = make_near_index(vec![NearInterval {
            start: 100,
            end: 110,
            group_id: None,
            strand: Strand::Plus,
        }]);
        let mut cfg = PrepareConfig::default();
        cfg.distance_sign = DistSign::Absolute;
        cfg.out_labels = vec!["win-direction".to_string()];

        let windows = vec![build_window("chr1", 50, 60, "A")]; // upstream

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result[0].label_tuples[0].near_side.as_deref(), Some("-"));
        Ok(())
    }

    #[test]
    fn far_hit_dropped_by_distance_max_absolute() {
        // Arrange: absolute distance beyond max drops window
        let near_index = make_near_index(vec![NearInterval {
            start: 0,
            end: 10,
            group_id: None,
            strand: Strand::Plus,
        }]);
        let mut cfg = PrepareConfig::default();
        cfg.distance_sign = DistSign::Absolute;
        cfg.distance_max = Some(5);
        let windows = vec![build_window("chr1", 100, 110, "A")]; // distance ~90

        // Act
        let result = apply_near_annotations(
            windows,
            &mut Some(near_index),
            &cfg,
            None,
            CoordinateSet::Resized,
        );

        // Assert
        assert!(result.is_empty());
    }

    #[test]
    fn left_edge_vs_right_edge_asymmetry() {
        // Arrange: window nearer right edge; left/right modes should differ
        let intervals = vec![NearInterval {
            start: 0,
            end: 100,
            group_id: None,
            strand: Strand::Plus,
        }];
        let mut chrom_left = NearChrom {
            intervals: intervals.clone(),
            cursor: 0,
        };
        let mut chrom_right = NearChrom {
            intervals,
            cursor: 0,
        };

        // Act
        let left = nearest_edge_distance(
            90,
            95,
            &mut chrom_left,
            &NearEdge::Left,
            &NearDirection::Both,
            true,
        )
        .expect("left hit");
        let right = nearest_edge_distance(
            90,
            95,
            &mut chrom_right,
            &NearEdge::Right,
            &NearDirection::Both,
            true,
        )
        .expect("right hit");

        // Assert: left edge distance is far, right edge is close
        if let NearestResult::Single(NearestDistance { distance, .. }) = left {
            assert_eq!(distance, 90);
        } else {
            panic!("left edge expected single hit");
        }
        if let NearestResult::Single(NearestDistance { distance, .. }) = right {
            assert_eq!(distance, -5);
        } else {
            panic!("right edge expected single hit");
        }
    }

    #[test]
    fn parse_distance_bins_rejects_reversed_range() {
        // Arrange / Act
        let result = parse_distance_bins(&vec!["bad:10-5".to_string()]);

        // Assert
        assert!(result.is_err());
    }

    #[test]
    fn parse_distance_bins_accepts_zero_length_range() -> Result<()> {
        // Arrange / Act
        let bins = parse_distance_bins(&vec!["zero:0-0".to_string()])?;

        // Assert
        assert_eq!(bins.match_label(0), Some("zero"));
        assert_eq!(bins.match_label(1), None);
        Ok(())
    }

    #[test]
    fn signed_bins_classify_upstream_overlap_and_downstream() -> Result<()> {
        // Arrange
        // Signed distances with bins that straddle zero; upstream negative, downstream positive, overlap zero
        let near_index = make_near_index(vec![NearInterval {
            start: 1000,
            end: 1010,
            group_id: None,
            strand: Strand::Plus,
        }]);
        let mut cfg = PrepareConfig::default();
        cfg.distance_sign = DistSign::Signed;
        cfg.distance_bins = Some(vec![
            "upstream:<0".to_string(),
            "at:0-0".to_string(),
            "downstream:>0".to_string(),
        ]);
        cfg.out_labels = vec!["win-direction".to_string(), "bin".to_string()];
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        // Upstream window at 800..820 -> nearest edge 1000 => distance -180
        let upstream = vec![build_window("chr1", 800, 820, "A")];
        let upstream_res = apply_near_annotations(
            upstream,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            upstream_res[0].label_tuples[0].bin.as_deref(),
            Some("upstream")
        );
        assert_eq!(
            upstream_res[0].label_tuples[0].near_side.as_deref(),
            Some("-")
        );

        // Overlap window 1005..1008 -> distance 0
        let overlap = vec![build_window("chr1", 1005, 1008, "A")];
        let overlap_res = apply_near_annotations(
            overlap,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(overlap_res[0].label_tuples[0].bin.as_deref(), Some("at"));
        assert_eq!(
            overlap_res[0].label_tuples[0].near_side.as_deref(),
            Some("=")
        );

        // Downstream window 1200..1210 -> distance +190
        let downstream = vec![build_window("chr1", 1200, 1210, "A")];
        let downstream_res = apply_near_annotations(
            downstream,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            downstream_res[0].label_tuples[0].bin.as_deref(),
            Some("downstream")
        );
        assert_eq!(
            downstream_res[0].label_tuples[0].near_side.as_deref(),
            Some("+")
        );
        Ok(())
    }

    #[test]
    fn signed_bins_with_negative_range_and_asymmetric_cutoffs() -> Result<()> {
        // Arrange
        // Bins cross zero with wider upstream span; expects:
        // dist1: upstream beyond -2500
        // prox: between -2500 and +500 (inclusive)
        // dist2: downstream beyond +500
        let near_index = make_near_index(vec![NearInterval {
            start: 5000,
            end: 5001,
            group_id: None,
            strand: Strand::Plus,
        }]);
        let mut cfg = PrepareConfig::default();
        cfg.distance_sign = DistSign::Signed;
        cfg.distance_bins = Some(vec![
            "prox:-2500-500".to_string(),
            "dist1:<-2500".to_string(),
            "dist2:>500".to_string(),
        ]);
        cfg.out_labels = vec!["win-direction".to_string(), "bin".to_string()];
        let bins = parse_distance_bins(cfg.distance_bins.as_ref().unwrap())?;

        // Very upstream: distance about -4000 should fall in dist1
        let upstream_far = vec![build_window("chr1", 900, 910, "A")];
        let upstream_far_res = apply_near_annotations(
            upstream_far,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            upstream_far_res[0].label_tuples[0].bin.as_deref(),
            Some("dist1")
        );
        assert_eq!(
            upstream_far_res[0].label_tuples[0].near_side.as_deref(),
            Some("-")
        );

        // Proximal window upstream but within -2500..500: distance about -1000 -> prox
        let prox_upstream = vec![build_window("chr1", 3800, 3810, "A")];
        let prox_upstream_res = apply_near_annotations(
            prox_upstream,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            prox_upstream_res[0].label_tuples[0].bin.as_deref(),
            Some("prox")
        );

        // Proximal downstream small: distance about +300 -> prox (still within +500)
        let prox_downstream = vec![build_window("chr1", 5300, 5310, "A")];
        let prox_downstream_res = apply_near_annotations(
            prox_downstream,
            &mut Some(near_index.clone()),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            prox_downstream_res[0].label_tuples[0].bin.as_deref(),
            Some("prox")
        );
        assert_eq!(
            prox_downstream_res[0].label_tuples[0].near_side.as_deref(),
            Some("+")
        );

        // Far downstream: distance about +1500 -> dist2
        let downstream_far = vec![build_window("chr1", 6500, 6510, "A")];
        let downstream_far_res = apply_near_annotations(
            downstream_far,
            &mut Some(near_index),
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );
        assert_eq!(
            downstream_far_res[0].label_tuples[0].bin.as_deref(),
            Some("dist2")
        );
        assert_eq!(
            downstream_far_res[0].label_tuples[0].near_side.as_deref(),
            Some("+")
        );
        Ok(())
    }
}
