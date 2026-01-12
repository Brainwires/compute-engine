// Unit tests for physics::gravitational_waves::detector
use computational_engine::compute::physics::gravitational_waves::detector::*;

use super::*;

    #[test]
    fn test_detector_creation() {
        let hanford = LIGODetector::hanford();
        let livingston = LIGODetector::livingston();

        assert_eq!(hanford.name, "LIGO Hanford");
        assert_eq!(livingston.name, "LIGO Livingston");
        assert!(hanford.latitude > 0.0); // Northern hemisphere
    }

    #[test]
    fn test_aligo_noise_curve() {
        let detector = LIGODetector::hanford();

        // Check noise at different frequencies
        let noise_low = detector.noise_asd(20.0);
        let noise_sweet = detector.noise_asd(100.0);
        let noise_high = detector.noise_asd(1000.0);

        // All should be finite and positive
        assert!(noise_low.is_finite() && noise_low > 0.0);
        assert!(noise_sweet.is_finite() && noise_sweet > 0.0);
        assert!(noise_high.is_finite() && noise_high > 0.0);
    }

    #[test]
    fn test_detector_network() {
        let network = DetectorNetwork::ligo_hl();

        assert_eq!(network.detectors.len(), 2);

        // Network sensitivity should be better than single detector
        let single_sens = network.detectors[0].noise_asd(100.0);
        let network_sens = network.network_sensitivity(100.0);

        // Network should be more sensitive (lower noise)
        assert!(network_sens < single_sens);
    }

    #[test]
    fn test_noise_models() {
        let mut detector = LIGODetector::hanford();

        detector.noise_model = NoiseModel::AdvancedLIGO;
        let aligo_noise = detector.noise_asd(100.0);

        detector.noise_model = NoiseModel::InitialLIGO;
        let iligo_noise = detector.noise_asd(100.0);

        detector.noise_model = NoiseModel::EinsteinTelescope;
        let et_noise = detector.noise_asd(100.0);

        // aLIGO should be better than iLIGO
        assert!(aligo_noise < iligo_noise);

        // ET should be best
        assert!(et_noise < aligo_noise);
    }

    #[test]
    fn test_antenna_pattern() {
        let detector = LIGODetector::hanford();

        let (f_plus, f_cross) = detector.antenna_pattern(0.0, 0.0, 0.0);

        // Antenna patterns should be between -1 and 1 (but simplified to 1)
        assert!(f_plus >= -1.0 && f_plus <= 1.0);
        assert!(f_cross >= -1.0 && f_cross <= 1.0);
    }
