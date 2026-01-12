use computational_engine::information_theory::*;

#[test]
fn test_shannon_entropy() {
    let data = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
    let result = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("shannon".to_string()),
        alpha: None,
    })
    .unwrap();

    assert!(result.entropy > 0.0);
    assert!(result.normalized_entropy <= 1.0);
    assert_eq!(result.base, 2.0);
}

#[test]
fn test_shannon_entropy_uniform() {
    // Uniform distribution should have maximum entropy
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let result = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("shannon".to_string()),
        alpha: None,
    })
    .unwrap();

    assert!(result.entropy > 0.0);
    assert!(result.normalized_entropy > 0.8); // Should be close to 1
}

#[test]
fn test_shannon_entropy_concentrated() {
    // Concentrated distribution should have low entropy
    let data = vec![1.0, 1.0, 1.0, 1.0, 1.0, 2.0];
    let result = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("shannon".to_string()),
        alpha: None,
    })
    .unwrap();

    assert!(result.entropy >= 0.0);
    assert!(result.normalized_entropy < 1.0);
}

#[test]
fn test_renyi_entropy() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let result = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("renyi".to_string()),
        alpha: Some(2.0),
    })
    .unwrap();

    assert!(result.entropy > 0.0);
    assert_eq!(result.entropy_type, "renyi");
}

#[test]
fn test_renyi_entropy_alpha_values() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result_alpha_05 = shannon_entropy(EntropyRequest {
        data: data.clone(),
        base: Some(2.0),
        entropy_type: Some("renyi".to_string()),
        alpha: Some(0.5),
    })
    .unwrap();

    let result_alpha_2 = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("renyi".to_string()),
        alpha: Some(2.0),
    })
    .unwrap();

    assert!(result_alpha_05.entropy > 0.0);
    assert!(result_alpha_2.entropy > 0.0);
}

#[test]
fn test_differential_entropy() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let result = shannon_entropy(EntropyRequest {
        data,
        base: Some(2.0),
        entropy_type: Some("differential".to_string()),
        alpha: None,
    })
    .unwrap();

    assert!(result.entropy > 0.0);
    assert_eq!(result.entropy_type, "differential");
}

#[test]
fn test_mutual_information() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0]; // y = 2x, should have high MI

    let result = mutual_information(MutualInfoRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    assert!(result.mutual_information > 0.0);
}

#[test]
fn test_mutual_information_independent() {
    let x = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
    let y = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];

    let result = mutual_information(MutualInfoRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    // Independent variables should have low MI
    assert!(result.mutual_information < 2.0);
}

#[test]
fn test_mutual_information_identical() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = x.clone();

    let result = mutual_information(MutualInfoRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    // Identical variables should have high MI
    assert!(result.mutual_information > 0.0);
}

#[test]
fn test_channel_capacity() {
    // Binary symmetric channel
    let transition_matrix = vec![vec![0.9, 0.1], vec![0.1, 0.9]];

    let result = channel_capacity(ChannelCapacityRequest {
        transition_matrix,
        input_distribution: None,
    })
    .unwrap();

    assert!(result.capacity > 0.0);
    assert!(result.capacity <= 1.0); // Capacity of binary channel is at most 1 bit
    assert_eq!(result.optimal_input_distribution.len(), 2);
}

#[test]
fn test_channel_capacity_noiseless() {
    // Noiseless channel (perfect transmission)
    let transition_matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    let result = channel_capacity(ChannelCapacityRequest {
        transition_matrix,
        input_distribution: None,
    })
    .unwrap();

    assert!((result.capacity - 1.0).abs() < 0.1); // Should be 1 bit for binary
}

#[test]
fn test_channel_capacity_custom_input() {
    let transition_matrix = vec![vec![0.8, 0.2], vec![0.3, 0.7]];

    let result = channel_capacity(ChannelCapacityRequest {
        transition_matrix,
        input_distribution: Some(vec![0.6, 0.4]),
    })
    .unwrap();

    assert!(result.capacity > 0.0);
}

#[test]
fn test_huffman_coding() {
    let symbols = vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ];
    let frequencies = vec![0.4, 0.3, 0.2, 0.1];

    let result = huffman_coding(HuffmanRequest {
        symbols,
        frequencies,
    })
    .unwrap();

    assert_eq!(result.codes.len(), 4);
    assert!(result.average_length > 0.0);
    assert!(result.efficiency > 0.0 && result.efficiency <= 1.0);
}

#[test]
fn test_huffman_coding_equal_frequencies() {
    let symbols = vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ];
    let frequencies = vec![0.25, 0.25, 0.25, 0.25];

    let result = huffman_coding(HuffmanRequest {
        symbols,
        frequencies,
    })
    .unwrap();

    assert_eq!(result.codes.len(), 4);
    // For equal frequencies, all codes should have similar length
    assert!(result.average_length > 1.5 && result.average_length < 2.5);
}

#[test]
fn test_huffman_coding_skewed() {
    let symbols = vec!["A".to_string(), "B".to_string()];
    let frequencies = vec![0.9, 0.1];

    let result = huffman_coding(HuffmanRequest {
        symbols,
        frequencies,
    })
    .unwrap();

    assert_eq!(result.codes.len(), 2);
    assert!(result.average_length < 2.0); // Should be efficient for skewed distribution
}

#[test]
fn test_kolmogorov_complexity() {
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

    let result = kolmogorov_complexity(KolmogorovRequest {
        data,
        compressor: Some("lz77".to_string()),
    })
    .unwrap();

    assert!(result.complexity_estimate > 0.0);
    assert!(result.compressed_size > 0);
    assert!(result.original_size > 0);
    assert!(result.compression_ratio > 0.0);
}

#[test]
fn test_kolmogorov_complexity_repetitive() {
    // Repetitive data should compress well
    let data = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

    let result = kolmogorov_complexity(KolmogorovRequest {
        data,
        compressor: Some("lz77".to_string()),
    })
    .unwrap();

    assert!(result.compression_ratio >= 1.0);
}

#[test]
fn test_kolmogorov_complexity_random() {
    // Random-looking data compresses poorly
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

    let result = kolmogorov_complexity(KolmogorovRequest {
        data,
        compressor: Some("lz77".to_string()),
    })
    .unwrap();

    assert!(result.compressed_size <= result.original_size);
}

#[test]
fn test_conditional_entropy() {
    let x = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
    let y = vec![1.0, 2.0, 2.0, 3.0, 3.0, 4.0];

    let result = conditional_entropy(ConditionalEntropyRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    assert!(result.conditional_entropy >= 0.0);
    assert!(result.normalized >= 0.0 && result.normalized <= 1.0);
}

#[test]
fn test_conditional_entropy_deterministic() {
    // If Y is fully determined by X, H(Y|X) should be low
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0]; // y = 2x

    let result = conditional_entropy(ConditionalEntropyRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    // Should be low since Y is determined by X
    assert!(result.conditional_entropy < 2.0);
}

#[test]
fn test_conditional_entropy_independent() {
    // If Y is independent of X, H(Y|X) = H(Y)
    let x = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
    let y = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];

    let result = conditional_entropy(ConditionalEntropyRequest {
        x,
        y,
        bins: Some(3),
    })
    .unwrap();

    assert!(result.conditional_entropy >= 0.0);
}

#[test]
fn test_relative_entropy() {
    let p = vec![0.5, 0.3, 0.2];
    let q = vec![0.4, 0.4, 0.2];

    let result = relative_entropy(RelativeEntropyRequest { p, q }).unwrap();

    assert!(result.kl_divergence >= 0.0);
    assert!(result.js_divergence >= 0.0);
}

#[test]
fn test_relative_entropy_identical() {
    let p = vec![0.25, 0.25, 0.25, 0.25];
    let q = vec![0.25, 0.25, 0.25, 0.25];

    let result = relative_entropy(RelativeEntropyRequest { p, q }).unwrap();

    assert!(result.kl_divergence.abs() < 1e-6); // Should be 0
    assert!(result.js_divergence.abs() < 1e-6); // Should be 0
}

#[test]
fn test_relative_entropy_symmetric() {
    let p = vec![0.7, 0.2, 0.1];
    let q = vec![0.3, 0.4, 0.3];

    let result1 = relative_entropy(RelativeEntropyRequest {
        p: p.clone(),
        q: q.clone(),
    })
    .unwrap();
    let result2 = relative_entropy(RelativeEntropyRequest { p: q, q: p }).unwrap();

    // KL divergence is not symmetric
    assert!((result1.kl_divergence - result2.kl_divergence).abs() > 0.001);
    // JS divergence is symmetric
    assert!((result1.js_divergence - result2.js_divergence).abs() < 1e-6);
}

#[test]
fn test_entropy_different_bases() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result_base2 = shannon_entropy(EntropyRequest {
        data: data.clone(),
        base: Some(2.0),
        entropy_type: Some("shannon".to_string()),
        alpha: None,
    })
    .unwrap();

    let result_base_e = shannon_entropy(EntropyRequest {
        data,
        base: Some(std::f64::consts::E),
        entropy_type: Some("shannon".to_string()),
        alpha: None,
    })
    .unwrap();

    // Different bases should give different entropy values
    assert!((result_base2.entropy - result_base_e.entropy).abs() > 0.1);
}

#[test]
fn test_huffman_compression_ratio() {
    let symbols = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let frequencies = vec![0.5, 0.3, 0.2];

    let result = huffman_coding(HuffmanRequest {
        symbols,
        frequencies,
    })
    .unwrap();

    // Compression ratio should be reasonable
    assert!(result.compression_ratio >= 1.0);
    assert!(result.compression_ratio <= 8.0);
}

#[test]
fn test_channel_capacity_three_inputs() {
    let transition_matrix = vec![
        vec![0.7, 0.2, 0.1],
        vec![0.1, 0.8, 0.1],
        vec![0.1, 0.1, 0.8],
    ];

    let result = channel_capacity(ChannelCapacityRequest {
        transition_matrix,
        input_distribution: None,
    })
    .unwrap();

    assert!(result.capacity > 0.0);
    assert_eq!(result.optimal_input_distribution.len(), 3);
}

#[test]
fn test_mutual_information_different_bins() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];

    let result_few_bins = mutual_information(MutualInfoRequest {
        x: x.clone(),
        y: y.clone(),
        bins: Some(2),
    })
    .unwrap();

    let result_many_bins = mutual_information(MutualInfoRequest {
        x,
        y,
        bins: Some(4),
    })
    .unwrap();

    // Both should be positive
    assert!(result_few_bins.mutual_information > 0.0);
    assert!(result_many_bins.mutual_information > 0.0);
}
