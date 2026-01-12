/**
 * Date/Time Mathematics Tests
 */
use computational_engine::datetime::*;

#[test]
fn test_add_days() {
    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: DateTimeParams {
            date1: Some("2025-01-01".to_string()),
            days: Some(30),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-31");
}

#[test]
fn test_add_months() {
    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: DateTimeParams {
            date1: Some("2025-01-31".to_string()),
            months: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    // Jan 31 + 1 month = Feb 28 (not 31)
    assert_eq!(result.value, "2025-02-28");
}

#[test]
fn test_add_years() {
    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: DateTimeParams {
            date1: Some("2024-02-29".to_string()), // Leap year
            years: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    // Feb 29, 2024 + 1 year = Feb 28, 2025 (not leap year)
    assert_eq!(result.value, "2025-02-28");
}

#[test]
fn test_date_difference() {
    let input = DateTimeInput {
        operation: DateTimeOperation::DateDifference,
        parameters: DateTimeParams {
            date1: Some("2025-01-01".to_string()),
            date2: Some("2025-01-31".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.numeric_value, Some(30.0));
    assert_eq!(result.unit, "days");
}

#[test]
fn test_subtract_days() {
    let input = DateTimeInput {
        operation: DateTimeOperation::SubtractInterval,
        parameters: DateTimeParams {
            date1: Some("2025-01-31".to_string()),
            days: Some(30),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-01");
}

#[test]
fn test_is_leap_year() {
    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: DateTimeParams {
            year: Some(2024),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "true");
    assert_eq!(result.numeric_value, Some(1.0));
}

#[test]
fn test_is_not_leap_year() {
    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: DateTimeParams {
            year: Some(2023),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "false");
}

#[test]
fn test_days_in_month_feb_leap() {
    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: DateTimeParams {
            year: Some(2024),
            month: Some(2),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.numeric_value, Some(29.0));
}

#[test]
fn test_days_in_month_feb_non_leap() {
    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: DateTimeParams {
            year: Some(2023),
            month: Some(2),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.numeric_value, Some(28.0));
}

#[test]
fn test_day_of_week() {
    let input = DateTimeInput {
        operation: DateTimeOperation::DayOfWeek,
        parameters: DateTimeParams {
            date1: Some("2025-10-05".to_string()), // Sunday
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "Sunday");
    assert_eq!(result.numeric_value, Some(0.0));
}

#[test]
fn test_business_days_simple() {
    let input = DateTimeInput {
        operation: DateTimeOperation::BusinessDays,
        parameters: DateTimeParams {
            date1: Some("2025-10-06".to_string()), // Monday
            date2: Some("2025-10-10".to_string()), // Friday
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    // Mon, Tue, Wed, Thu, Fri = 5 business days
    assert_eq!(result.numeric_value, Some(5.0));
}

#[test]
fn test_age_calculation() {
    let input = DateTimeInput {
        operation: DateTimeOperation::AgeAtDate,
        parameters: DateTimeParams {
            date1: Some("2000-01-01".to_string()), // Birthdate
            date2: Some("2025-01-01".to_string()), // Target date
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.numeric_value, Some(25.0)); // Approximately 25 years
}

#[test]
fn test_week_number() {
    let input = DateTimeInput {
        operation: DateTimeOperation::WeekNumber,
        parameters: DateTimeParams {
            date1: Some("2025-01-15".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    assert!(result.numeric_value.unwrap() >= 1.0 && result.numeric_value.unwrap() <= 53.0);
}

#[test]
fn test_add_mixed_intervals() {
    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: DateTimeParams {
            date1: Some("2025-01-01".to_string()),
            days: Some(15),
            months: Some(2),
            years: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_datetime(input).unwrap();
    // 2025-01-01 + 1 year = 2026-01-01
    // 2026-01-01 + 2 months = 2026-03-01
    // 2026-03-01 + 15 days = 2026-03-16
    assert_eq!(result.value, "2026-03-16");
}

#[test]
fn test_century_leap_year() {
    // 2000 was a leap year (divisible by 400)
    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: DateTimeParams {
            year: Some(2000),
            ..Default::default()
        },
    };
    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "true");

    // 1900 was NOT a leap year (divisible by 100 but not 400)
    let input2 = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: DateTimeParams {
            year: Some(1900),
            ..Default::default()
        },
    };
    let result2 = calculate_datetime(input2).unwrap();
    assert_eq!(result2.value, "false");
}
