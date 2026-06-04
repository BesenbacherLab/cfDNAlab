mod tests_transitions_frequency_calculations {
    use crate::commands::transitions::transitions::compute_transition_frequencies;
    use anyhow::Result;
    use ndarray::{array, s};

    #[test]
    fn computes_first_order_frequencies() -> Result<()> {
        // Arrange
        let counts = array![[[2.0, 1.0, 4.0, 3.0]]];
        let motifs = vec![
            "AA".to_string(),
            "AC".to_string(),
            "CA".to_string(),
            "CC".to_string(),
        ];

        // Act
        let freqs = compute_transition_frequencies(&counts, 1, &motifs)?;

        // Assert
        assert!((freqs[[0, 0, 0]] - (2.0 / 3.0)).abs() < 1e-9);
        assert!((freqs[[0, 0, 1]] - (1.0 / 3.0)).abs() < 1e-9);
        assert!((freqs[[0, 0, 2]] - (4.0 / 7.0)).abs() < 1e-9);
        assert!((freqs[[0, 0, 3]] - (3.0 / 7.0)).abs() < 1e-9);
        Ok(())
    }

    #[test]
    fn computes_second_order_frequencies_with_multiple_positions() -> Result<()> {
        // Arrange
        let counts = array![
            [[4.0, 6.0, 2.0, 0.0], [0.0, 0.0, 5.0, 5.0],],
            [[1.0, 1.0, 0.0, 4.0], [0.0, 0.0, 0.0, 0.0],],
        ];
        let motifs = vec![
            "AAT".to_string(),
            "AAC".to_string(),
            "CAT".to_string(),
            "CAC".to_string(),
        ];

        // Act
        let freqs = compute_transition_frequencies(&counts, 2, &motifs)?;

        // Assert
        assert!((freqs[[0, 0, 0]] - 0.4).abs() < 1e-9);
        assert!((freqs[[0, 0, 1]] - 0.6).abs() < 1e-9);
        assert!((freqs[[0, 0, 2]] - 1.0).abs() < 1e-9);
        assert!((freqs[[0, 0, 3]] - 0.0).abs() < 1e-9);

        assert!((freqs[[0, 1, 2]] - 0.5).abs() < 1e-9);
        assert!((freqs[[0, 1, 3]] - 0.5).abs() < 1e-9);

        assert!((freqs[[1, 0, 0]] - 0.5).abs() < 1e-9);
        assert!((freqs[[1, 0, 1]] - 0.5).abs() < 1e-9);
        assert!((freqs[[1, 0, 2]] - 0.0).abs() < 1e-9);
        assert!((freqs[[1, 0, 3]] - 1.0).abs() < 1e-9);

        assert!(freqs.slice(s![1, 1, ..]).iter().all(|v| v.abs() < 1e-9));
        Ok(())
    }

    #[test]
    fn errors_when_motif_axis_mismatches_counts() {
        // Arrange
        let counts = array![[[1.0, 2.0]]];
        let motifs = vec!["AA".to_string()];

        // Act
        let result = compute_transition_frequencies(&counts, 1, &motifs);

        // Assert
        assert!(result.is_err());
    }
}
