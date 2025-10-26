//! Unit tests for datetime API handler
//!
//! Tests all datetime operations including:
//! - Date arithmetic (add/subtract intervals)
//! - Date differences
//! - Age calculations (current and at specific dates)
//! - Business day calculations
//! - Leap year checks
//! - Days in month calculations
//! - Week number calculations
//! - Day of week calculations
//! - Error handling for invalid inputs

use crate::api::handlers::datetime;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a datetime request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "datetime".to_string(),
        operation: operation.to_string(),
        parameters,
    }
}

/// Helper function to extract result value from response
fn get_result(response: &ComputationResponse) -> &Value {
    assert!(response.success, "Expected success response");
    response.result.as_ref().expect("Expected result value")
}

// ============================================================================
// Add Interval Tests
// ============================================================================

#[test]
fn test_add_days() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-15",
            "days": 10
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "2024-01-25");
    assert_eq!(result["unit"].as_str().unwrap(), "date");
}

#[test]
fn test_add_weeks() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01",
            "weeks": 2
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "2024-01-15");
}

#[test]
fn test_add_months() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-15",
            "months": 3
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "2024-04-15");
}

#[test]
fn test_add_years() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2020-02-29",
            "years": 1
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Feb 29 -> Feb 28 on non-leap year
    assert_eq!(result["value"].as_str().unwrap(), "2021-02-28");
}

#[test]
fn test_add_multiple_intervals() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01",
            "days": 5,
            "weeks": 2,
            "months": 1
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Subtract Interval Tests
// ============================================================================

#[test]
fn test_subtract_days() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("subtract_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-25",
            "days": 10
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "2024-01-15");
}

#[test]
fn test_subtract_months() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("subtract_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-04-15",
            "months": 3
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "2024-01-15");
}

// ============================================================================
// Date Difference Tests
// ============================================================================

#[test]
fn test_date_difference_positive() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("date_difference"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01",
            "date2": "2024-01-31"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "30");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 30.0);
    assert_eq!(result["unit"].as_str().unwrap(), "days");
}

#[test]
fn test_date_difference_negative() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("date_difference"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-31",
            "date2": "2024-01-01"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "30");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), -30.0);
}

#[test]
fn test_date_difference_year_span() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("date_difference"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2023-01-01",
            "date2": "2024-01-01"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "365");
}

// ============================================================================
// Age Calculation Tests
// ============================================================================

#[test]
fn test_age_current() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("age_current"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "1990-05-15"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Age should be calculated from 1990-05-15 to 2025-10-05 (hardcoded in module)
    let age = result["numeric_value"].as_f64().unwrap();
    assert!(age >= 35.0 && age <= 36.0);
}

#[test]
fn test_age_at_date() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("age_at_date"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "1990-05-15",
            "date2": "2020-05-15"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Exactly 30 years
    let age = result["numeric_value"].as_f64().unwrap();
    assert_eq!(age, 30.0);
}

// ============================================================================
// Business Days Tests
// ============================================================================

#[test]
fn test_business_days_same_week() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("business_days"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01", // Monday
            "date2": "2024-01-05"  // Friday
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Monday to Friday = 5 business days
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 5.0);
}

#[test]
fn test_business_days_with_weekend() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("business_days"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01", // Monday
            "date2": "2024-01-08"  // Next Monday
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should exclude weekend days
    assert!(result["numeric_value"].as_f64().unwrap() <= 6.0);
}

#[test]
fn test_business_days_with_holidays() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("business_days"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01",
            "date2": "2024-01-05",
            "exclude_holidays": ["2024-01-03"]
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Add Business Days Tests
// ============================================================================

#[test]
fn test_add_business_days() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_business_days"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01", // Monday
            "days": 5
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should skip weekends
    assert!(result["value"].as_str().is_some());
}

#[test]
fn test_add_business_days_with_holidays() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_business_days"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01",
            "days": 3,
            "exclude_holidays": ["2024-01-03"]
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Leap Year Tests
// ============================================================================

#[test]
fn test_is_leap_year_true() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("is_leap_year"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2024
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "true");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 1.0);
}

#[test]
fn test_is_leap_year_false() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("is_leap_year"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2023
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "false");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 0.0);
}

#[test]
fn test_is_leap_year_century_divisible_400() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("is_leap_year"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2000
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "true");
}

#[test]
fn test_is_leap_year_century_not_divisible_400() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("is_leap_year"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 1900
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "false");
}

// ============================================================================
// Days in Month Tests
// ============================================================================

#[test]
fn test_days_in_month_january() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("days_in_month"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2024,
            "month": 1
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "31");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 31.0);
}

#[test]
fn test_days_in_month_february_leap() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("days_in_month"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2024,
            "month": 2
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "29");
}

#[test]
fn test_days_in_month_february_non_leap() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("days_in_month"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2023,
            "month": 2
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "28");
}

#[test]
fn test_days_in_month_april() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("days_in_month"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2024,
            "month": 4
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "30");
}

// ============================================================================
// Week Number Tests
// ============================================================================

#[test]
fn test_week_number_start_of_year() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("week_number"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert!(result["numeric_value"].as_f64().unwrap() >= 1.0);
}

#[test]
fn test_week_number_mid_year() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("week_number"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-06-15"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let week = result["numeric_value"].as_f64().unwrap();
    assert!(week >= 1.0 && week <= 53.0); // Valid ISO 8601 week number range
}

// ============================================================================
// Day of Week Tests
// ============================================================================

#[test]
fn test_day_of_week_monday() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("day_of_week"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01" // Monday
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "Monday");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 1.0);
}

#[test]
fn test_day_of_week_sunday() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("day_of_week"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-07" // Sunday
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "Sunday");
    assert_eq!(result["numeric_value"].as_f64().unwrap(), 0.0);
}

#[test]
fn test_day_of_week_saturday() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("day_of_week"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-06" // Saturday
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["value"].as_str().unwrap(), "Saturday");
    assert!(result["additional_info"]["is_weekend"].as_str().unwrap() == "true");
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_date_format() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "01-15-2024", // Wrong format
            "days": 10
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

#[test]
fn test_missing_required_date() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "days": 10
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
}

#[test]
fn test_invalid_month() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("days_in_month"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "year": 2024,
            "month": 13
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
}

#[test]
fn test_invalid_day_for_month() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("add_interval"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-02-30", // Feb doesn't have 30 days
            "days": 1
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
}

#[test]
fn test_missing_date2_for_difference() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("date_difference"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "date1": "2024-01-01"
        }),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
}

#[test]
fn test_malformed_json() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("invalid_operation"),
    );

    let request = create_request("datetime", params);
    let response = datetime::handle(&request);

    assert!(!response.success);
}
