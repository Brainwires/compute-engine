// Unit tests for datetime::mod
use computational_engine::datetime::mod::*;

use super::*;

    // Basic Date Tests (5 existing tests to keep)
    #[test]
    fn test_date_parsing() {
        let date = Date::from_iso_string("2025-10-05").unwrap();
        assert_eq!(date.year, 2025);
        assert_eq!(date.month, 10);
        assert_eq!(date.day, 5);
    }

    #[test]
    fn test_leap_year() {
        assert!(Date::is_leap_year(2024));
        assert!(!Date::is_leap_year(2023));
        assert!(Date::is_leap_year(2000));
        assert!(!Date::is_leap_year(1900));
    }

    #[test]
    fn test_add_days() {
        let date = Date::new(2025, 1, 1).unwrap();
        let new_date = date.add_days(31);
        assert_eq!(new_date.month, 2);
        assert_eq!(new_date.day, 1);
    }

    #[test]
    fn test_date_difference() {
        let date1 = Date::new(2025, 1, 1).unwrap();
        let date2 = Date::new(2025, 1, 31).unwrap();
        assert_eq!(date1.days_between(&date2), 30);
    }

    #[test]
    fn test_day_of_week() {
        let date = Date::new(2025, 10, 5).unwrap(); // This is a Sunday
        assert_eq!(date.day_of_week(), 0); // 0 = Sunday
    }

    // Add Interval Tests (3 tests)
    #[test]
    fn test_add_interval_days() {
        let params = DateTimeParams {
            date1: Some("2025-01-01".to_string()),
            days: Some(10),
            ..Default::default()
        };
        let result = add_interval(&params).unwrap();
        assert_eq!(result.value, "2025-01-11");
    }

    #[test]
    fn test_add_interval_months() {
        let params = DateTimeParams {
            date1: Some("2025-01-31".to_string()),
            months: Some(1),
            ..Default::default()
        };
        let result = add_interval(&params).unwrap();
        // Jan 31 + 1 month = Feb 28 (day capped to valid)
        assert_eq!(result.value, "2025-02-28");
    }

    #[test]
    fn test_add_interval_years() {
        let params = DateTimeParams {
            date1: Some("2024-02-29".to_string()), // Leap year
            years: Some(1),
            ..Default::default()
        };
        let result = add_interval(&params).unwrap();
        // Feb 29 2024 + 1 year = Feb 28 2025 (non-leap)
        assert_eq!(result.value, "2025-02-28");
    }

    // Subtract Interval Tests (2 tests)
    #[test]
    fn test_subtract_interval_days() {
        let params = DateTimeParams {
            date1: Some("2025-01-15".to_string()),
            days: Some(10),
            ..Default::default()
        };
        let result = subtract_interval(&params).unwrap();
        assert_eq!(result.value, "2025-01-05");
    }

    #[test]
    fn test_subtract_interval_months() {
        let params = DateTimeParams {
            date1: Some("2025-03-31".to_string()),
            months: Some(1),
            ..Default::default()
        };
        let result = subtract_interval(&params).unwrap();
        // Mar 31 - 1 month = Feb 28
        assert_eq!(result.value, "2025-02-28");
    }

    // Date Difference Tests (2 tests)
    #[test]
    fn test_date_difference_positive() {
        let params = DateTimeParams {
            date1: Some("2025-01-01".to_string()),
            date2: Some("2025-01-31".to_string()),
            ..Default::default()
        };
        let result = date_difference(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 30.0);
    }

    #[test]
    fn test_date_difference_negative() {
        let params = DateTimeParams {
            date1: Some("2025-01-31".to_string()),
            date2: Some("2025-01-01".to_string()),
            ..Default::default()
        };
        let result = date_difference(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), -30.0);
    }

    // Age Tests (2 tests)
    #[test]
    fn test_age_at_date() {
        let params = DateTimeParams {
            date1: Some("2000-01-01".to_string()), // Birthdate
            date2: Some("2025-01-01".to_string()), // Target date
            ..Default::default()
        };
        let result = age_at_date(&params).unwrap();
        // 25 years (9131 days / 365)
        assert!(result.numeric_value.unwrap() >= 24.0 && result.numeric_value.unwrap() <= 26.0);
    }

    #[test]
    fn test_age_current() {
        let params = DateTimeParams {
            date1: Some("1990-10-05".to_string()), // Birthdate
            ..Default::default()
        };
        let result = age_current(&params).unwrap();
        // Should be around 35 years old (as of 2025-10-05)
        assert!(result.numeric_value.unwrap() >= 34.0 && result.numeric_value.unwrap() <= 36.0);
    }

    // Business Days Tests (3 tests)
    #[test]
    fn test_business_days_one_week() {
        let params = DateTimeParams {
            date1: Some("2025-10-06".to_string()), // Monday
            date2: Some("2025-10-10".to_string()), // Friday
            ..Default::default()
        };
        let result = business_days(&params).unwrap();
        // Mon-Fri = 5 business days
        assert_eq!(result.numeric_value.unwrap(), 5.0);
    }

    #[test]
    fn test_business_days_with_weekend() {
        let params = DateTimeParams {
            date1: Some("2025-10-03".to_string()), // Friday
            date2: Some("2025-10-06".to_string()), // Monday
            ..Default::default()
        };
        let result = business_days(&params).unwrap();
        // Fri, Sat, Sun, Mon = 2 business days (Fri + Mon)
        assert_eq!(result.numeric_value.unwrap(), 2.0);
    }

    #[test]
    fn test_business_days_same_day() {
        let params = DateTimeParams {
            date1: Some("2025-10-06".to_string()),
            date2: Some("2025-10-06".to_string()),
            ..Default::default()
        };
        let result = business_days(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 1.0); // Same day counts as 1
    }

    // Is Leap Year Tests (2 tests)
    #[test]
    fn test_is_leap_year_operation() {
        let params = DateTimeParams {
            year: Some(2024),
            ..Default::default()
        };
        let result = is_leap_year(&params).unwrap();
        assert_eq!(result.value, "true");
        assert_eq!(result.numeric_value.unwrap(), 1.0);
    }

    #[test]
    fn test_is_not_leap_year() {
        let params = DateTimeParams {
            year: Some(2023),
            ..Default::default()
        };
        let result = is_leap_year(&params).unwrap();
        assert_eq!(result.value, "false");
        assert_eq!(result.numeric_value.unwrap(), 0.0);
    }

    // Days in Month Tests (3 tests)
    #[test]
    fn test_days_in_month_february_leap() {
        let params = DateTimeParams {
            year: Some(2024),
            month: Some(2),
            ..Default::default()
        };
        let result = days_in_month(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 29.0);
    }

    #[test]
    fn test_days_in_month_february_non_leap() {
        let params = DateTimeParams {
            year: Some(2023),
            month: Some(2),
            ..Default::default()
        };
        let result = days_in_month(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 28.0);
    }

    #[test]
    fn test_days_in_month_31_days() {
        let params = DateTimeParams {
            year: Some(2025),
            month: Some(1), // January
            ..Default::default()
        };
        let result = days_in_month(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 31.0);
    }

    // Week Number Tests (2 tests)
    #[test]
    fn test_week_number_first_week() {
        let params = DateTimeParams {
            date1: Some("2025-01-05".to_string()),
            ..Default::default()
        };
        let result = week_number(&params).unwrap();
        assert!(result.numeric_value.unwrap() >= 1.0);
    }

    #[test]
    fn test_week_number_middle_of_year() {
        let params = DateTimeParams {
            date1: Some("2025-07-01".to_string()),
            ..Default::default()
        };
        let result = week_number(&params).unwrap();
        // Week should be in first half of year (1-53)
        assert!(result.numeric_value.unwrap() >= 1.0 && result.numeric_value.unwrap() <= 53.0);
    }

    // Day of Week Tests (3 tests)
    #[test]
    fn test_day_of_week_operation_monday() {
        let params = DateTimeParams {
            date1: Some("2025-10-06".to_string()), // Monday
            ..Default::default()
        };
        let result = day_of_week(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 1.0); // 1 = Monday
        assert!(result.value.contains("Monday"));
    }

    #[test]
    fn test_day_of_week_sunday() {
        let params = DateTimeParams {
            date1: Some("2025-10-05".to_string()), // Sunday
            ..Default::default()
        };
        let result = day_of_week(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 0.0); // 0 = Sunday
        assert!(result.value.contains("Sunday"));
    }

    #[test]
    fn test_day_of_week_weekend() {
        let params = DateTimeParams {
            date1: Some("2025-10-04".to_string()), // Saturday
            ..Default::default()
        };
        let result = day_of_week(&params).unwrap();
        assert_eq!(result.numeric_value.unwrap(), 6.0); // 6 = Saturday
        let info = result.additional_info.unwrap();
        assert_eq!(info.get("is_weekend").unwrap(), "true");
    }

    // Add Business Days Tests (2 tests)
    #[test]
    fn test_add_business_days_within_week() {
        let params = DateTimeParams {
            date1: Some("2025-10-06".to_string()), // Monday
            days: Some(3),
            ..Default::default()
        };
        let result = add_business_days(&params).unwrap();
        // Mon + 3 business days = Thursday
        assert_eq!(result.value, "2025-10-09");
    }

    #[test]
    fn test_add_business_days_crossing_weekend() {
        let params = DateTimeParams {
            date1: Some("2025-10-03".to_string()), // Friday
            days: Some(3),
            ..Default::default()
        };
        let result = add_business_days(&params).unwrap();
        // Fri + 3 business days = Wed (skip Sat, Sun, Mon, Tue, Wed)
        assert_eq!(result.value, "2025-10-08");
    }
