use super::*;

// Test AddInterval operation
#[test]
fn test_add_days_to_date() {
    let params = DateTimeParams {
        date1: Some("2025-01-15".to_string()),
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-25");
    assert_eq!(result.numeric_value, Some(10.0));
}

#[test]
fn test_add_weeks_to_date() {
    let params = DateTimeParams {
        date1: Some("2025-01-01".to_string()),
        weeks: Some(2),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-15");
    assert_eq!(result.numeric_value, Some(14.0));
}

#[test]
fn test_add_months_to_date() {
    let params = DateTimeParams {
        date1: Some("2025-01-31".to_string()),
        months: Some(1),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    // January 31 + 1 month should give February 28 (not a leap year)
    assert_eq!(result.value, "2025-02-28");
}

#[test]
fn test_add_years_to_date() {
    let params = DateTimeParams {
        date1: Some("2020-02-29".to_string()),
        years: Some(1),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    // Feb 29 2020 + 1 year = Feb 28 2021 (not a leap year)
    assert_eq!(result.value, "2021-02-28");
}

#[test]
fn test_add_combined_intervals() {
    let params = DateTimeParams {
        date1: Some("2025-01-01".to_string()),
        days: Some(5),
        months: Some(2),
        years: Some(1),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    // 2025-01-01 + 1 year = 2026-01-01, + 2 months = 2026-03-01, + 5 days = 2026-03-06
    assert_eq!(result.value, "2026-03-06");
}

// Test SubtractInterval operation
#[test]
fn test_subtract_days_from_date() {
    let params = DateTimeParams {
        date1: Some("2025-01-15".to_string()),
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::SubtractInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-05");
}

#[test]
fn test_subtract_months_from_date() {
    let params = DateTimeParams {
        date1: Some("2025-03-31".to_string()),
        months: Some(1),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::SubtractInterval,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    // March 31 - 1 month = February 28 (not a leap year)
    assert_eq!(result.value, "2025-02-28");
}

// Test DateDifference operation
#[test]
fn test_date_difference_forward() {
    let params = DateTimeParams {
        date1: Some("2025-01-01".to_string()),
        date2: Some("2025-01-31".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DateDifference,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "30");
    assert_eq!(result.numeric_value, Some(30.0));
    assert_eq!(result.unit, "days");
}

#[test]
fn test_date_difference_backward() {
    let params = DateTimeParams {
        date1: Some("2025-01-31".to_string()),
        date2: Some("2025-01-01".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DateDifference,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "30");
    assert_eq!(result.numeric_value, Some(-30.0));
}

// Test IsLeapYear operation
#[test]
fn test_is_leap_year_true() {
    let params = DateTimeParams {
        year: Some(2024),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "true");
    assert_eq!(result.numeric_value, Some(1.0));
}

#[test]
fn test_is_leap_year_false() {
    let params = DateTimeParams {
        year: Some(2025),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "false");
    assert_eq!(result.numeric_value, Some(0.0));
}

#[test]
fn test_is_leap_year_century_not_divisible_by_400() {
    let params = DateTimeParams {
        year: Some(1900),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "false");
}

#[test]
fn test_is_leap_year_century_divisible_by_400() {
    let params = DateTimeParams {
        year: Some(2000),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::IsLeapYear,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "true");
}

// Test DaysInMonth operation
#[test]
fn test_days_in_month_january() {
    let params = DateTimeParams {
        year: Some(2025),
        month: Some(1),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "31");
    assert_eq!(result.numeric_value, Some(31.0));
}

#[test]
fn test_days_in_month_february_non_leap() {
    let params = DateTimeParams {
        year: Some(2025),
        month: Some(2),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "28");
}

#[test]
fn test_days_in_month_february_leap() {
    let params = DateTimeParams {
        year: Some(2024),
        month: Some(2),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "29");
}

#[test]
fn test_days_in_month_april() {
    let params = DateTimeParams {
        year: Some(2025),
        month: Some(4),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DaysInMonth,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "30");
}

// Test DayOfWeek operation
#[test]
fn test_day_of_week_known_date() {
    // January 1, 2025 is a Wednesday
    let params = DateTimeParams {
        date1: Some("2025-01-01".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DayOfWeek,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "Wednesday");
    assert_eq!(result.numeric_value, Some(3.0));
}

#[test]
fn test_day_of_week_sunday() {
    // January 5, 2025 is a Sunday
    let params = DateTimeParams {
        date1: Some("2025-01-05".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DayOfWeek,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "Sunday");
    assert_eq!(result.numeric_value, Some(0.0));
}

#[test]
fn test_day_of_week_saturday() {
    // January 4, 2025 is a Saturday
    let params = DateTimeParams {
        date1: Some("2025-01-04".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::DayOfWeek,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "Saturday");
    assert_eq!(result.numeric_value, Some(6.0));
}

// Test BusinessDays operation
#[test]
fn test_business_days_no_weekends() {
    // Monday to Friday (5 business days)
    let params = DateTimeParams {
        date1: Some("2025-01-06".to_string()), // Monday
        date2: Some("2025-01-10".to_string()), // Friday
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::BusinessDays,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "5");
}

#[test]
fn test_business_days_with_weekends() {
    // Monday to Monday (5 business days, excluding weekend)
    let params = DateTimeParams {
        date1: Some("2025-01-06".to_string()), // Monday
        date2: Some("2025-01-13".to_string()), // Monday next week
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::BusinessDays,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "6");
}

// Test AddBusinessDays operation
#[test]
fn test_add_business_days_no_weekend() {
    // Starting Monday, add 3 business days = Thursday
    let params = DateTimeParams {
        date1: Some("2025-01-06".to_string()), // Monday
        days: Some(3),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddBusinessDays,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-09"); // Thursday
}

#[test]
fn test_add_business_days_skip_weekend() {
    // Starting Friday, add 3 business days = Wednesday (skip weekend)
    let params = DateTimeParams {
        date1: Some("2025-01-10".to_string()), // Friday
        days: Some(3),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddBusinessDays,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert_eq!(result.value, "2025-01-15"); // Wednesday
}

// Test WeekNumber operation
#[test]
fn test_week_number_january() {
    let params = DateTimeParams {
        date1: Some("2025-01-15".to_string()),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::WeekNumber,
        parameters: params,
    };

    let result = calculate_datetime(input).unwrap();
    assert!(result.numeric_value.unwrap() >= 1.0);
    assert!(result.numeric_value.unwrap() <= 53.0);
}

// Test error handling
#[test]
fn test_invalid_date_format() {
    let params = DateTimeParams {
        date1: Some("2025/01/15".to_string()), // Wrong format
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input);
    assert!(result.is_err());
}

#[test]
fn test_invalid_date_values() {
    let params = DateTimeParams {
        date1: Some("2025-13-01".to_string()), // Invalid month
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input);
    assert!(result.is_err());
}

#[test]
fn test_invalid_day_for_month() {
    let params = DateTimeParams {
        date1: Some("2025-02-30".to_string()), // February doesn't have 30 days
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input);
    assert!(result.is_err());
}

#[test]
fn test_missing_required_parameter() {
    let params = DateTimeParams {
        // Missing date1
        days: Some(10),
        ..Default::default()
    };

    let input = DateTimeInput {
        operation: DateTimeOperation::AddInterval,
        parameters: params,
    };

    let result = calculate_datetime(input);
    assert!(result.is_err());
}
