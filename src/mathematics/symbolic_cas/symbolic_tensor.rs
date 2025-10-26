//! Symbolic tensor operations
//!
//! Support for arbitrary-rank tensors with symbolic expressions

use super::{Expr, SymbolicError, SymbolicMatrix, SymbolicResult};
use std::fmt;

/// Index type for tensor notation
#[derive(Debug, Clone, PartialEq)]
pub enum IndexType {
    /// Contravariant index (upper): T^i
    Contravariant,
    /// Covariant index (lower): T_i
    Covariant,
}

/// A symbolic tensor of arbitrary rank
#[derive(Debug, Clone, PartialEq)]
pub struct SymbolicTensor {
    /// The rank of the tensor (number of indices)
    rank: usize,

    /// Dimension of each index
    dimensions: Vec<usize>,

    /// Index types (contravariant or covariant)
    index_types: Vec<IndexType>,

    /// Flat storage of tensor components as symbolic expressions
    /// For rank-2: [T_00, T_01, ..., T_10, T_11, ...]
    data: Vec<Expr>,
}

impl SymbolicTensor {
    /// Create a new tensor from data and index types
    pub fn new(
        dimensions: Vec<usize>,
        index_types: Vec<IndexType>,
        data: Vec<Expr>,
    ) -> SymbolicResult<Self> {
        let rank = dimensions.len();

        if rank != index_types.len() {
            return Err(SymbolicError::InvalidOperation(format!(
                "Rank mismatch: {} dimensions but {} index types",
                rank,
                index_types.len()
            )));
        }

        let expected_size: usize = dimensions.iter().product();
        if data.len() != expected_size {
            return Err(SymbolicError::InvalidOperation(format!(
                "Data size mismatch: expected {}, got {}",
                expected_size,
                data.len()
            )));
        }

        Ok(SymbolicTensor {
            rank,
            dimensions,
            index_types,
            data,
        })
    }

    /// Create a scalar (rank-0 tensor)
    pub fn scalar(value: Expr) -> Self {
        SymbolicTensor {
            rank: 0,
            dimensions: vec![],
            index_types: vec![],
            data: vec![value],
        }
    }

    /// Create a vector (rank-1 tensor)
    pub fn vector(components: Vec<Expr>, index_type: IndexType) -> SymbolicResult<Self> {
        let n = components.len();
        Ok(SymbolicTensor {
            rank: 1,
            dimensions: vec![n],
            index_types: vec![index_type],
            data: components,
        })
    }

    /// Create a rank-2 tensor from a symbolic matrix
    pub fn from_matrix(matrix: &SymbolicMatrix, index_types: [IndexType; 2]) -> Self {
        let rows = matrix.rows();
        let cols = matrix.cols();
        let mut data = Vec::with_capacity(rows * cols);

        for i in 0..rows {
            for j in 0..cols {
                data.push(matrix.get(i, j).unwrap().clone());
            }
        }

        SymbolicTensor {
            rank: 2,
            dimensions: vec![rows, cols],
            index_types: index_types.to_vec(),
            data,
        }
    }

    /// Convert rank-2 tensor to matrix (if possible)
    pub fn to_matrix(&self) -> SymbolicResult<SymbolicMatrix> {
        if self.rank != 2 {
            return Err(SymbolicError::InvalidOperation(format!(
                "Cannot convert rank-{} tensor to matrix",
                self.rank
            )));
        }

        let rows = self.dimensions[0];
        let cols = self.dimensions[1];
        let mut matrix_data = vec![vec![Expr::num(0); cols]; rows];

        for i in 0..rows {
            for j in 0..cols {
                let idx = i * cols + j;
                matrix_data[i][j] = self.data[idx].clone();
            }
        }

        SymbolicMatrix::new(matrix_data)
    }

    /// Get tensor rank
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Get dimensions
    pub fn dimensions(&self) -> &[usize] {
        &self.dimensions
    }

    /// Get index types
    pub fn index_types(&self) -> &[IndexType] {
        &self.index_types
    }

    /// Get tensor data (components)
    pub fn data(&self) -> &[Expr] {
        &self.data
    }

    /// Get component by multi-index
    pub fn get(&self, indices: &[usize]) -> SymbolicResult<&Expr> {
        if indices.len() != self.rank {
            return Err(SymbolicError::InvalidOperation(format!(
                "Expected {} indices, got {}",
                self.rank,
                indices.len()
            )));
        }

        for (i, &idx) in indices.iter().enumerate() {
            if idx >= self.dimensions[i] {
                return Err(SymbolicError::InvalidOperation(format!(
                    "Index {} out of bounds for dimension {}",
                    idx, i
                )));
            }
        }

        let flat_index = self.compute_flat_index(indices);
        Ok(&self.data[flat_index])
    }

    /// Set component by multi-index
    pub fn set(&mut self, indices: &[usize], value: Expr) -> SymbolicResult<()> {
        if indices.len() != self.rank {
            return Err(SymbolicError::InvalidOperation(format!(
                "Expected {} indices, got {}",
                self.rank,
                indices.len()
            )));
        }

        for (i, &idx) in indices.iter().enumerate() {
            if idx >= self.dimensions[i] {
                return Err(SymbolicError::InvalidOperation(format!(
                    "Index {} out of bounds for dimension {}",
                    idx, i
                )));
            }
        }

        let flat_index = self.compute_flat_index(indices);
        self.data[flat_index] = value;
        Ok(())
    }

    /// Compute flat array index from multi-index
    fn compute_flat_index(&self, indices: &[usize]) -> usize {
        let mut flat_idx = 0;
        let mut stride = 1;

        for i in (0..self.rank).rev() {
            flat_idx += indices[i] * stride;
            stride *= self.dimensions[i];
        }

        flat_idx
    }

    /// Tensor contraction: sum over paired indices
    /// For rank-2 tensor, this is the trace
    /// For higher ranks, this reduces the rank by 2
    pub fn contract(&self, index1: usize, index2: usize) -> SymbolicResult<Self> {
        if index1 >= self.rank || index2 >= self.rank {
            return Err(SymbolicError::InvalidOperation(
                "Index out of bounds".to_string(),
            ));
        }

        if index1 == index2 {
            return Err(SymbolicError::InvalidOperation(
                "Cannot contract the same index".to_string(),
            ));
        }

        if self.dimensions[index1] != self.dimensions[index2] {
            return Err(SymbolicError::InvalidOperation(
                "Contracted indices must have the same dimension".to_string(),
            ));
        }

        // Check that one is contravariant and one is covariant
        if self.index_types[index1] == self.index_types[index2] {
            return Err(SymbolicError::InvalidOperation(
                "Contracted indices must be of opposite types (one upper, one lower)".to_string(),
            ));
        }

        // For rank-2, contraction gives a scalar (trace)
        if self.rank == 2 {
            let mut sum = Expr::num(0);
            let n = self.dimensions[0];
            for i in 0..n {
                let idx = i * self.dimensions[1] + i;
                sum = Expr::add(sum, self.data[idx].clone());
            }
            return Ok(SymbolicTensor::scalar(sum));
        }

        // For higher ranks, we need to build a new tensor with rank-2
        let (idx1, idx2) = if index1 < index2 {
            (index1, index2)
        } else {
            (index2, index1)
        };

        // Build new dimensions and index types (remove contracted indices)
        let mut new_dimensions = Vec::new();
        let mut new_index_types = Vec::new();
        for i in 0..self.rank {
            if i != idx1 && i != idx2 {
                new_dimensions.push(self.dimensions[i]);
                new_index_types.push(self.index_types[i].clone());
            }
        }

        let new_rank = self.rank - 2;
        if new_rank == 0 {
            // Contraction results in a scalar
            return self.contract_to_scalar(idx1, idx2);
        }

        // Calculate size of new tensor
        let new_size: usize = new_dimensions.iter().product();
        let mut new_data = vec![Expr::num(0); new_size];

        // The contraction dimension
        let contract_dim = self.dimensions[idx1];

        // Build new tensor by summing over contracted indices
        self.compute_contraction(&mut new_data, &new_dimensions, idx1, idx2, contract_dim)?;

        Ok(SymbolicTensor {
            rank: new_rank,
            dimensions: new_dimensions,
            index_types: new_index_types,
            data: new_data,
        })
    }

    /// Helper: contract to scalar (for rank > 2)
    fn contract_to_scalar(&self, idx1: usize, idx2: usize) -> SymbolicResult<Self> {
        let contract_dim = self.dimensions[idx1];
        let mut sum = Expr::num(0);

        // Iterate over all possible values of the contracted index
        for k in 0..contract_dim {
            // Build multi-index with both contracted indices set to k
            let mut indices = vec![0; self.rank];
            indices[idx1] = k;
            indices[idx2] = k;

            // For other indices, they don't exist in result, so we sum over all
            // This handles the general case where we're summing T^i_i
            let flat_idx = self.compute_flat_index(&indices);
            sum = Expr::add(sum, self.data[flat_idx].clone());
        }

        Ok(SymbolicTensor::scalar(sum))
    }

    /// Helper: compute contraction for higher-rank tensors
    fn compute_contraction(
        &self,
        new_data: &mut [Expr],
        new_dimensions: &[usize],
        idx1: usize,
        idx2: usize,
        contract_dim: usize,
    ) -> SymbolicResult<()> {
        let new_rank = new_dimensions.len();

        // Iterate over all elements of the new tensor
        let new_size = new_data.len();
        for flat_new_idx in 0..new_size {
            // Convert flat index to multi-index for new tensor
            let mut new_multi_idx = vec![0; new_rank];
            let mut temp = flat_new_idx;
            for i in (0..new_rank).rev() {
                new_multi_idx[i] = temp % new_dimensions[i];
                temp /= new_dimensions[i];
            }

            // Sum over the contracted index
            let mut sum = Expr::num(0);
            for k in 0..contract_dim {
                // Map new multi-index back to old multi-index
                let mut old_multi_idx = Vec::with_capacity(self.rank);
                let mut j = 0;
                for i in 0..self.rank {
                    if i == idx1 || i == idx2 {
                        old_multi_idx.push(k);
                    } else {
                        old_multi_idx.push(new_multi_idx[j]);
                        j += 1;
                    }
                }

                let old_flat_idx = self.compute_flat_index(&old_multi_idx);
                sum = Expr::add(sum, self.data[old_flat_idx].clone());
            }

            new_data[flat_new_idx] = sum;
        }

        Ok(())
    }

    /// Outer product of two tensors
    pub fn outer_product(&self, other: &SymbolicTensor) -> SymbolicResult<Self> {
        let new_rank = self.rank + other.rank;
        let mut new_dimensions = self.dimensions.clone();
        new_dimensions.extend_from_slice(&other.dimensions);

        let mut new_index_types = self.index_types.clone();
        new_index_types.extend_from_slice(&other.index_types);

        let new_size = self.data.len() * other.data.len();
        let mut new_data = Vec::with_capacity(new_size);

        for elem1 in &self.data {
            for elem2 in &other.data {
                new_data.push(Expr::mul(elem1.clone(), elem2.clone()));
            }
        }

        Ok(SymbolicTensor {
            rank: new_rank,
            dimensions: new_dimensions,
            index_types: new_index_types,
            data: new_data,
        })
    }

    /// Raise/lower an index using a metric tensor
    pub fn raise_index(&self, metric: &SymbolicMatrix, index: usize) -> SymbolicResult<Self> {
        if self.rank != 1 {
            return Err(SymbolicError::InvalidOperation(
                "Index raising only implemented for vectors".to_string(),
            ));
        }

        if self.index_types[index] != IndexType::Covariant {
            return Err(SymbolicError::InvalidOperation(
                "Can only raise covariant indices".to_string(),
            ));
        }

        let n = self.dimensions[0];
        if metric.rows() != n || metric.cols() != n {
            return Err(SymbolicError::InvalidOperation(
                "Metric dimension mismatch".to_string(),
            ));
        }

        // v^i = g^ij v_j (sum over j)
        let mut new_components = Vec::with_capacity(n);
        for i in 0..n {
            let mut sum = Expr::num(0);
            for j in 0..n {
                let g_ij = metric.get(i, j).unwrap().clone();
                let v_j = self.data[j].clone();
                sum = Expr::add(sum, Expr::mul(g_ij, v_j));
            }
            new_components.push(sum);
        }

        SymbolicTensor::vector(new_components, IndexType::Contravariant)
    }

    /// Lower an index using a metric tensor
    pub fn lower_index(&self, metric: &SymbolicMatrix, index: usize) -> SymbolicResult<Self> {
        if self.rank != 1 {
            return Err(SymbolicError::InvalidOperation(
                "Index lowering only implemented for vectors".to_string(),
            ));
        }

        if self.index_types[index] != IndexType::Contravariant {
            return Err(SymbolicError::InvalidOperation(
                "Can only lower contravariant indices".to_string(),
            ));
        }

        let n = self.dimensions[0];
        if metric.rows() != n || metric.cols() != n {
            return Err(SymbolicError::InvalidOperation(
                "Metric dimension mismatch".to_string(),
            ));
        }

        // v_i = g_ij v^j (sum over j)
        let mut new_components = Vec::with_capacity(n);
        for i in 0..n {
            let mut sum = Expr::num(0);
            for j in 0..n {
                let g_ij = metric.get(i, j).unwrap().clone();
                let v_j = self.data[j].clone();
                sum = Expr::add(sum, Expr::mul(g_ij, v_j));
            }
            new_components.push(sum);
        }

        SymbolicTensor::vector(new_components, IndexType::Covariant)
    }
}

impl fmt::Display for SymbolicTensor {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.rank {
            0 => write!(f, "Scalar: {}", self.data[0]),
            1 => {
                let index_symbol = match self.index_types[0] {
                    IndexType::Contravariant => "^i",
                    IndexType::Covariant => "_i",
                };
                writeln!(f, "Vector T{}:", index_symbol)?;
                for (i, component) in self.data.iter().enumerate() {
                    writeln!(f, "  [{}] = {}", i, component)?;
                }
                Ok(())
            }
            2 => {
                let idx1 = match self.index_types[0] {
                    IndexType::Contravariant => "^i",
                    IndexType::Covariant => "_i",
                };
                let idx2 = match self.index_types[1] {
                    IndexType::Contravariant => "^j",
                    IndexType::Covariant => "_j",
                };
                writeln!(f, "Rank-2 Tensor T{}{}:", idx1, idx2)?;
                let rows = self.dimensions[0];
                let cols = self.dimensions[1];
                for i in 0..rows {
                    write!(f, "  [")?;
                    for j in 0..cols {
                        if j > 0 {
                            write!(f, ", ")?;
                        }
                        let idx = i * cols + j;
                        write!(f, "{}", self.data[idx])?;
                    }
                    writeln!(f, "]")?;
                }
                Ok(())
            }
            _ => {
                writeln!(f, "Rank-{} Tensor:", self.rank)?;
                writeln!(f, "Dimensions: {:?}", self.dimensions)?;
                writeln!(f, "Components: {} total", self.data.len())?;
                Ok(())
            }
        }
    }
}

