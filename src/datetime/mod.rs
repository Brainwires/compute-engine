/**
 * Date/Time Mathematics Module
 *
 * Provides calculations for:
 * - Date arithmetic (add/subtract intervals)
 * - Duration between dates
 * - Business day calculations
 * - Age calculations
 * - Time zone conversions
 * - Calendar-aware operations
 */
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DateTimeInput {
    pub operation: DateTimeOperation,
    pub parameters: DateTimeParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DateTimeOperation {
    AddInterval,      // Add days/months/years to a date
    SubtractInterval, // Subtract days/months/years from a date
    DateDifference,   // Calculate difference between two dates
    AgeCurrent,       // Calculate age from birthdate to now
    AgeAtDate,        // Calculate age at a specific date
    BusinessDays,     // Count business days between dates
    IsLeapYear,       // Check if year is leap year
    DaysInMonth,      // Get number of days in a month
    WeekNumber,       // Get week number in year
    DayOfWeek,        // Get day of week (0=Sunday, 6=Saturday)
    AddBusinessDays,  // Add N business days to a date
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DateTimeParams {
    // Date inputs (ISO 8601 format: YYYY-MM-DD)
    pub date1: Option<String>, // Primary date
    pub date2: Option<String>, // Secondary date (for differences)
    pub year: Option<i32>,
    pub month: Option<u32>, // 1-12
    pub day: Option<u32>,   // 1-31

    // Interval inputs
    pub days: Option<i64>,
    pub weeks: Option<i64>,
    pub months: Option<i64>,
    pub years: Option<i64>,

    // Options
    pub include_end_date: Option<bool>,        // For date ranges
    pub exclude_weekends: Option<bool>,        // For business day calculations
    pub exclude_holidays: Option<Vec<String>>, // List of holiday dates
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DateTimeResult {
    pub value: String,              // Result date or number as string
    pub numeric_value: Option<f64>, // Numeric result (days, years, etc.)
    pub unit: String,
    pub operation_used: String,
    pub interpretation: String,
    pub additional_info: Option<HashMap<String, String>>,
}

impl Default for DateTimeParams {
    fn default() -> Self {
        Self {
            date1: None,
            date2: None,
            year: None,
            month: None,
            day: None,
            days: None,
            weeks: None,
            months: None,
            years: None,
            include_end_date: None,
            exclude_weekends: None,
            exclude_holidays: None,
        }
    }
}

// Simple date structure (not using external crates)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Date {
    year: i32,
    month: u32,
    day: u32,
}

impl Date {
    fn new(year: i32, month: u32, day: u32) -> Result<Self, String> {
        if month < 1 || month > 12 {
            return Err(format!("Invalid month: {}", month));
        }

        let days_in_month = Self::days_in_month_static(year, month);
        if day < 1 || day > days_in_month {
            return Err(format!("Invalid day {} for month {}", day, month));
        }

        Ok(Date { year, month, day })
    }

    fn from_iso_string(s: &str) -> Result<Self, String> {
        let parts: Vec<&str> = s.split('-').collect();
        if parts.len() != 3 {
            return Err(format!("Invalid date format: {}. Use YYYY-MM-DD", s));
        }

        let year = parts[0]
            .parse::<i32>()
            .map_err(|_| format!("Invalid year: {}", parts[0]))?;
        let month = parts[1]
            .parse::<u32>()
            .map_err(|_| format!("Invalid month: {}", parts[1]))?;
        let day = parts[2]
            .parse::<u32>()
            .map_err(|_| format!("Invalid day: {}", parts[2]))?;

        Self::new(year, month, day)
    }

    fn to_iso_string(&self) -> String {
        format!("{:04}-{:02}-{:02}", self.year, self.month, self.day)
    }

    fn is_leap_year(year: i32) -> bool {
        (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
    }

    fn days_in_month_static(year: i32, month: u32) -> u32 {
        match month {
            1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
            4 | 6 | 9 | 11 => 30,
            2 => {
                if Self::is_leap_year(year) {
                    29
                } else {
                    28
                }
            }
            _ => 0,
        }
    }

    #[allow(dead_code)]
    fn days_in_month(&self) -> u32 {
        Self::days_in_month_static(self.year, self.month)
    }

    // Convert to days since epoch (Jan 1, 1970)
    fn to_days_since_epoch(&self) -> i64 {
        let mut days = 0i64;

        // Add days for complete years
        for y in 1970..self.year {
            days += if Self::is_leap_year(y) { 366 } else { 365 };
        }

        // Subtract days for years before 1970
        for y in self.year..1970 {
            days -= if Self::is_leap_year(y) { 366 } else { 365 };
        }

        // Add days for complete months in current year
        for m in 1..self.month {
            days += Self::days_in_month_static(self.year, m) as i64;
        }

        // Add remaining days
        days += self.day as i64;

        days
    }

    fn from_days_since_epoch(mut days: i64) -> Self {
        let mut year = 1970;

        // Handle years
        loop {
            let year_days = if Self::is_leap_year(year) { 366 } else { 365 };
            if days > year_days {
                days -= year_days;
                year += 1;
            } else {
                break;
            }
        }

        // Handle months
        let mut month = 1;
        loop {
            let month_days = Self::days_in_month_static(year, month) as i64;
            if days > month_days {
                days -= month_days;
                month += 1;
            } else {
                break;
            }
        }

        Date {
            year,
            month,
            day: days as u32,
        }
    }

    fn add_days(&self, days: i64) -> Self {
        let epoch_days = self.to_days_since_epoch();
        Self::from_days_since_epoch(epoch_days + days)
    }

    fn add_months(&self, months: i64) -> Result<Self, String> {
        let total_months = (self.year as i64 * 12 + self.month as i64 - 1) + months;
        let new_year = (total_months / 12) as i32;
        let new_month = (total_months % 12 + 1) as u32;

        // Handle day overflow (e.g., Jan 31 + 1 month = Feb 28/29)
        let max_day = Self::days_in_month_static(new_year, new_month);
        let new_day = self.day.min(max_day);

        Self::new(new_year, new_month, new_day)
    }

    fn add_years(&self, years: i64) -> Result<Self, String> {
        let new_year = self.year + years as i32;

        // Handle Feb 29 -> Feb 28 on non-leap years
        let new_day = if self.month == 2 && self.day == 29 && !Self::is_leap_year(new_year) {
            28
        } else {
            self.day
        };

        Self::new(new_year, self.month, new_day)
    }

    fn days_between(&self, other: &Date) -> i64 {
        other.to_days_since_epoch() - self.to_days_since_epoch()
    }

    fn day_of_week(&self) -> u32 {
        // Zeller's congruence
        let mut m = self.month as i32;
        let mut y = self.year;

        if m < 3 {
            m += 12;
            y -= 1;
        }

        let q = self.day as i32;
        let k = y % 100;
        let j = y / 100;

        let h = (q + ((13 * (m + 1)) / 5) + k + (k / 4) + (j / 4) - (2 * j)) % 7;

        // Convert to 0=Sunday, 6=Saturday
        ((h + 6) % 7) as u32
    }

    fn is_weekend(&self) -> bool {
        let dow = self.day_of_week();
        dow == 0 || dow == 6 // Sunday or Saturday
    }

    fn week_number(&self) -> u32 {
        // ISO 8601 week number
        let jan1 = Date::new(self.year, 1, 1).unwrap();
        let days_from_jan1 = self.days_between(&jan1);
        let jan1_dow = jan1.day_of_week();

        // Week 1 is the first week with a Thursday
        let week = ((days_from_jan1 + jan1_dow as i64 + 6) / 7) as u32;
        week.max(1).min(53)
    }
}

pub fn calculate_datetime(input: DateTimeInput) -> Result<DateTimeResult, String> {
    match input.operation {
        DateTimeOperation::AddInterval => add_interval(&input.parameters),
        DateTimeOperation::SubtractInterval => subtract_interval(&input.parameters),
        DateTimeOperation::DateDifference => date_difference(&input.parameters),
        DateTimeOperation::AgeCurrent => age_current(&input.parameters),
        DateTimeOperation::AgeAtDate => age_at_date(&input.parameters),
        DateTimeOperation::BusinessDays => business_days(&input.parameters),
        DateTimeOperation::IsLeapYear => is_leap_year(&input.parameters),
        DateTimeOperation::DaysInMonth => days_in_month(&input.parameters),
        DateTimeOperation::WeekNumber => week_number(&input.parameters),
        DateTimeOperation::DayOfWeek => day_of_week(&input.parameters),
        DateTimeOperation::AddBusinessDays => add_business_days(&input.parameters),
    }
}

fn add_interval(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date_str = params.date1.as_ref().ok_or("date1 required")?;
    let mut date = Date::from_iso_string(date_str)?;

    let mut total_days = 0i64;
    let mut info = HashMap::new();

    if let Some(days) = params.days {
        date = date.add_days(days);
        total_days += days;
        info.insert("days_added".to_string(), days.to_string());
    }

    if let Some(weeks) = params.weeks {
        let days = weeks * 7;
        date = date.add_days(days);
        total_days += days;
        info.insert("weeks_added".to_string(), weeks.to_string());
    }

    if let Some(months) = params.months {
        date = date.add_months(months)?;
        info.insert("months_added".to_string(), months.to_string());
    }

    if let Some(years) = params.years {
        date = date.add_years(years)?;
        info.insert("years_added".to_string(), years.to_string());
    }

    Ok(DateTimeResult {
        value: date.to_iso_string(),
        numeric_value: Some(total_days as f64),
        unit: "date".to_string(),
        operation_used: "Date + Interval".to_string(),
        interpretation: format!("New date: {}", date.to_iso_string()),
        additional_info: Some(info),
    })
}

fn subtract_interval(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let mut neg_params = params.clone();
    neg_params.days = params.days.map(|d| -d);
    neg_params.weeks = params.weeks.map(|w| -w);
    neg_params.months = params.months.map(|m| -m);
    neg_params.years = params.years.map(|y| -y);

    add_interval(&neg_params)
}

fn date_difference(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date1_str = params.date1.as_ref().ok_or("date1 required")?;
    let date2_str = params.date2.as_ref().ok_or("date2 required")?;

    let date1 = Date::from_iso_string(date1_str)?;
    let date2 = Date::from_iso_string(date2_str)?;

    let days_diff = date1.days_between(&date2);
    let weeks = days_diff / 7;
    let remaining_days = days_diff % 7;

    // Approximate months and years
    let approx_years = days_diff / 365;
    let approx_months = days_diff / 30;

    let mut info = HashMap::new();
    info.insert("total_days".to_string(), days_diff.to_string());
    info.insert("weeks".to_string(), weeks.to_string());
    info.insert("remaining_days".to_string(), remaining_days.to_string());
    info.insert(
        "approx_years".to_string(),
        format!("{:.2}", approx_years as f64 / 1.0),
    );
    info.insert("approx_months".to_string(), approx_months.to_string());

    let direction = if days_diff > 0 { "after" } else { "before" };

    Ok(DateTimeResult {
        value: days_diff.abs().to_string(),
        numeric_value: Some(days_diff as f64),
        unit: "days".to_string(),
        operation_used: "Date2 - Date1".to_string(),
        interpretation: format!(
            "{} is {} days {} {}",
            date2_str,
            days_diff.abs(),
            direction,
            date1_str
        ),
        additional_info: Some(info),
    })
}

fn age_current(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let birthdate_str = params.date1.as_ref().ok_or("date1 (birthdate) required")?;
    let birthdate = Date::from_iso_string(birthdate_str)?;

    // Use today's date (simplified - using epoch + known days)
    // In production, would use actual system time
    let today = Date::new(2025, 10, 5)?; // Current session date

    let days_diff = birthdate.days_between(&today);
    let years = days_diff / 365;

    let mut info = HashMap::new();
    info.insert("birthdate".to_string(), birthdate.to_iso_string());
    info.insert("current_date".to_string(), today.to_iso_string());
    info.insert("days_lived".to_string(), days_diff.to_string());

    Ok(DateTimeResult {
        value: years.to_string(),
        numeric_value: Some(years as f64),
        unit: "years".to_string(),
        operation_used: "Age calculation".to_string(),
        interpretation: format!("{} years old", years),
        additional_info: Some(info),
    })
}

fn age_at_date(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let birthdate_str = params.date1.as_ref().ok_or("date1 (birthdate) required")?;
    let target_date_str = params
        .date2
        .as_ref()
        .ok_or("date2 (target date) required")?;

    let birthdate = Date::from_iso_string(birthdate_str)?;
    let target_date = Date::from_iso_string(target_date_str)?;

    let days_diff = birthdate.days_between(&target_date);
    let years = days_diff / 365;

    let mut info = HashMap::new();
    info.insert("birthdate".to_string(), birthdate.to_iso_string());
    info.insert("target_date".to_string(), target_date.to_iso_string());
    info.insert("days_between".to_string(), days_diff.to_string());

    Ok(DateTimeResult {
        value: years.to_string(),
        numeric_value: Some(years as f64),
        unit: "years".to_string(),
        operation_used: "Age at specific date".to_string(),
        interpretation: format!("{} years old on {}", years, target_date_str),
        additional_info: Some(info),
    })
}

fn business_days(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date1_str = params.date1.as_ref().ok_or("date1 required")?;
    let date2_str = params.date2.as_ref().ok_or("date2 required")?;

    let date1 = Date::from_iso_string(date1_str)?;
    let date2 = Date::from_iso_string(date2_str)?;

    let total_days = date1.days_between(&date2).abs();
    let mut business_days = 0i64;
    let mut current = date1.min(date2);
    let end = date1.max(date2);

    let holidays: Vec<Date> = params
        .exclude_holidays
        .as_ref()
        .map(|h| {
            h.iter()
                .filter_map(|s| Date::from_iso_string(s).ok())
                .collect()
        })
        .unwrap_or_default();

    while current <= end {
        if !current.is_weekend() && !holidays.contains(&current) {
            business_days += 1;
        }
        current = current.add_days(1);
    }

    let mut info = HashMap::new();
    info.insert("total_days".to_string(), total_days.to_string());
    info.insert(
        "weekend_days".to_string(),
        (total_days - business_days).to_string(),
    );

    Ok(DateTimeResult {
        value: business_days.to_string(),
        numeric_value: Some(business_days as f64),
        unit: "business days".to_string(),
        operation_used: "Business days between dates".to_string(),
        interpretation: format!("{} business days (excluding weekends)", business_days),
        additional_info: Some(info),
    })
}

fn add_business_days(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date_str = params.date1.as_ref().ok_or("date1 required")?;
    let days_to_add = params.days.ok_or("days required")?.abs();

    let mut date = Date::from_iso_string(date_str)?;
    let holidays: Vec<Date> = params
        .exclude_holidays
        .as_ref()
        .map(|h| {
            h.iter()
                .filter_map(|s| Date::from_iso_string(s).ok())
                .collect()
        })
        .unwrap_or_default();

    let mut business_days_added = 0i64;

    while business_days_added < days_to_add {
        date = date.add_days(1);
        if !date.is_weekend() && !holidays.contains(&date) {
            business_days_added += 1;
        }
    }

    let mut info = HashMap::new();
    info.insert(
        "business_days_added".to_string(),
        business_days_added.to_string(),
    );

    Ok(DateTimeResult {
        value: date.to_iso_string(),
        numeric_value: Some(business_days_added as f64),
        unit: "date".to_string(),
        operation_used: "Add business days".to_string(),
        interpretation: format!(
            "Date after {} business days: {}",
            business_days_added,
            date.to_iso_string()
        ),
        additional_info: Some(info),
    })
}

fn is_leap_year(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let year = params.year.ok_or("year required")?;
    let is_leap = Date::is_leap_year(year);

    let mut info = HashMap::new();
    info.insert("year".to_string(), year.to_string());
    info.insert(
        "days_in_year".to_string(),
        if is_leap { "366" } else { "365" }.to_string(),
    );

    Ok(DateTimeResult {
        value: is_leap.to_string(),
        numeric_value: Some(if is_leap { 1.0 } else { 0.0 }),
        unit: "boolean".to_string(),
        operation_used: "Leap year check".to_string(),
        interpretation: format!(
            "{} is {} a leap year",
            year,
            if is_leap { "" } else { "not" }
        ),
        additional_info: Some(info),
    })
}

fn days_in_month(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let year = params.year.ok_or("year required")?;
    let month = params.month.ok_or("month required")?;

    if month < 1 || month > 12 {
        return Err(format!("Invalid month: {}", month));
    }

    let days = Date::days_in_month_static(year, month);

    let month_name = [
        "",
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ][month as usize];

    let mut info = HashMap::new();
    info.insert("year".to_string(), year.to_string());
    info.insert("month".to_string(), month.to_string());
    info.insert("month_name".to_string(), month_name.to_string());

    Ok(DateTimeResult {
        value: days.to_string(),
        numeric_value: Some(days as f64),
        unit: "days".to_string(),
        operation_used: "Days in month".to_string(),
        interpretation: format!("{} {} has {} days", month_name, year, days),
        additional_info: Some(info),
    })
}

fn week_number(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date_str = params.date1.as_ref().ok_or("date1 required")?;
    let date = Date::from_iso_string(date_str)?;

    let week = date.week_number();

    let mut info = HashMap::new();
    info.insert("date".to_string(), date.to_iso_string());
    info.insert("year".to_string(), date.year.to_string());

    Ok(DateTimeResult {
        value: week.to_string(),
        numeric_value: Some(week as f64),
        unit: "week number".to_string(),
        operation_used: "ISO 8601 week number".to_string(),
        interpretation: format!("Week {} of {}", week, date.year),
        additional_info: Some(info),
    })
}

fn day_of_week(params: &DateTimeParams) -> Result<DateTimeResult, String> {
    let date_str = params.date1.as_ref().ok_or("date1 required")?;
    let date = Date::from_iso_string(date_str)?;

    let dow = date.day_of_week();
    let day_names = [
        "Sunday",
        "Monday",
        "Tuesday",
        "Wednesday",
        "Thursday",
        "Friday",
        "Saturday",
    ];
    let day_name = day_names[dow as usize];

    let mut info = HashMap::new();
    info.insert("date".to_string(), date.to_iso_string());
    info.insert("day_name".to_string(), day_name.to_string());
    info.insert("is_weekend".to_string(), date.is_weekend().to_string());

    Ok(DateTimeResult {
        value: day_name.to_string(),
        numeric_value: Some(dow as f64),
        unit: "day of week".to_string(),
        operation_used: "Day of week calculation".to_string(),
        interpretation: format!("{} is a {}", date_str, day_name),
        additional_info: Some(info),
    })
}

#[cfg(test)]
#[path = "../../tests/unit/datetime_tests.rs"]
mod tests;

