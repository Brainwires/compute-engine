//! Information theory operations
//!
//! Provides entropy calculations, mutual information, channel capacity, and related operations

use crate::engine::*;
use crate::specialized::information_theory::*;

/// Compute information theory operations
pub fn compute_information(op: &InformationOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    let result_json = match op {
        InformationOp::Entropy => {
            let req: EntropyRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse entropy request: {}", e))?;

            let result =
                shannon_entropy(req).map_err(|e| format!("Entropy calculation error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::MutualInfo => {
            let req: MutualInfoRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse mutual info request: {}", e))?;

            let result =
                mutual_information(req).map_err(|e| format!("Mutual information error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::ChannelCapacity => {
            let req: ChannelCapacityRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse channel capacity request: {}", e))?;

            let result =
                channel_capacity(req).map_err(|e| format!("Failed to serialize result: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::Huffman => {
            let req: HuffmanRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse Huffman request: {}", e))?;

            let result = huffman_coding(req).map_err(|e| format!("Huffman coding error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::Kolmogorov => {
            let req: KolmogorovRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse Kolmogorov request: {}", e))?;

            let result = kolmogorov_complexity(req)
                .map_err(|e| format!("Kolmogorov complexity error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::ConditionalEntropy => {
            let req: ConditionalEntropyRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse conditional entropy request: {}", e))?;

            let result = conditional_entropy(req)
                .map_err(|e| format!("Conditional entropy error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        InformationOp::KLDivergence => {
            let req: RelativeEntropyRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse KL divergence request: {}", e))?;

            let result =
                relative_entropy(req).map_err(|e| format!("KL divergence error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}
