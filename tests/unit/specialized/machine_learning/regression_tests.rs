// Unit tests for specialized::machine_learning::regression
use computational_engine::specialized::machine_learning::*;

use super::*;

    #[test]
    fn test_linear_regression_simple() {
        // y = 2x (perfect linear relationship)
        let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![2.0, 4.0, 6.0, 8.0];

        let model = linear_regression(&x, &y).unwrap();

        // Coefficient should be close to 2.0
        assert!((model.coefficients[0] - 2.0).abs() < 0.01);

        // Intercept should be close to 0.0
        assert!(model.intercept.abs() < 0.01);

        // R² should be close to 1.0
        assert!(model.r_squared.unwrap() > 0.99);
    }

    #[test]
    fn test_linear_regression_with_intercept() {
        // y = 3x + 5
        let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![8.0, 11.0, 14.0, 17.0];

        let model = linear_regression(&x, &y).unwrap();

        assert!((model.coefficients[0] - 3.0).abs() < 0.01);
        assert!((model.intercept - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_linear_regression_predict() {
        let x = vec![vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![2.0, 4.0, 6.0];

        let model = linear_regression(&x, &y).unwrap();

        let x_test = vec![vec![5.0], vec![6.0]];
        let predictions = linear_regression_predict(&model, &x_test);

        // Should predict approximately [10.0, 12.0]
        assert!((predictions[0] - 10.0).abs() < 0.1);
        assert!((predictions[1] - 12.0).abs() < 0.1);
    }

    #[test]
    fn test_ridge_regression() {
        let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![2.0, 4.0, 6.0, 8.0];

        let model = ridge_regression(&x, &y, 1.0).unwrap();

        // With regularization, coefficients should be slightly smaller
        assert!(model.coefficients[0] < 2.0);
        assert!(model.coefficients[0] > 1.5);
    }

    #[test]
    fn test_logistic_regression_simple() {
        // Linearly separable data
        let x = vec![
            vec![1.0],
            vec![2.0],
            vec![3.0],
            vec![4.0],
            vec![5.0],
            vec![6.0],
        ];
        let y = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

        let model = logistic_regression(&x, &y, 0.1, 1000, 1e-6).unwrap();

        // Should achieve high accuracy
        assert!(model.accuracy.unwrap() > 0.8);

        // Coefficient should be positive (higher x → higher probability)
        assert!(model.coefficients[0] > 0.0);
    }

    #[test]
    fn test_logistic_regression_predict() {
        let x = vec![vec![1.0], vec![5.0]];
        let y = vec![0.0, 1.0];

        let model = logistic_regression(&x, &y, 0.1, 1000, 1e-6).unwrap();

        let x_test = vec![vec![0.5], vec![6.0]];
        let predictions = logistic_regression_predict(&model, &x_test);

        // Low x should predict 0, high x should predict 1
        assert_eq!(predictions[0], 0.0);
        assert_eq!(predictions[1], 1.0);
    }

    #[test]
    fn test_sigmoid() {
        assert!((sigmoid(0.0) - 0.5).abs() < 1e-10);
        assert!(sigmoid(100.0) > 0.99);
        assert!(sigmoid(-100.0) < 0.01);
    }

    #[test]
    fn test_r_squared() {
        let y_true = vec![1.0, 2.0, 3.0, 4.0];
        let y_pred = vec![1.1, 1.9, 3.2, 3.8];

        let r2 = compute_r_squared(&y_true, &y_pred);

        // Should be close to 1.0 for good fit
        assert!(r2 > 0.95);
    }

    #[test]
    fn test_mean_squared_error() {
        let y_true = vec![1.0, 2.0, 3.0];
        let y_pred = vec![1.0, 2.0, 3.0];

        let mse = mean_squared_error(&y_true, &y_pred);
        assert_eq!(mse, 0.0);
    }

    #[test]
    fn test_mean_absolute_error() {
        let y_true = vec![1.0, 2.0, 3.0];
        let y_pred = vec![1.5, 2.5, 3.5];

        let mae = mean_absolute_error(&y_true, &y_pred);
        assert_eq!(mae, 0.5);
    }
