use clap::{Parser, Subcommand};
use computational_engine::*;
use serde_json::json;

#[derive(Parser)]
#[command(name = "brainwires-compute-engine")]
#[command(about = "Unified computational engine for mathematics and physics", long_about = None)]
#[command(version)]
#[command(disable_help_subcommand = true)]
#[command(arg_required_else_help = false)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Display help for tools and their operations (e.g., "help compute", "help compute physics")
    Help {
        /// Tool name (solve, differentiate, integrate, analyze, simulate, compute, transform, fieldtheory, sample, optimize)
        /// If not specified, shows all operations across all tools
        tool: Option<String>,
        /// Category within the tool (e.g., "physics", "matrix", "tensor")
        category: Option<String>,
    },
    /// Tensor calculus and Einstein equations
    Tensor {
        #[command(subcommand)]
        operation: TensorOps,
    },
    /// Advanced calculus operations
    Calculus {
        #[command(subcommand)]
        operation: CalculusOps,
    },
    /// Fluid dynamics simulations
    Fluid {
        #[command(subcommand)]
        operation: FluidOps,
    },
    /// Signal processing
    Signal {
        #[command(subcommand)]
        operation: SignalOps,
    },
    /// Stochastic processes
    Stochastic {
        #[command(subcommand)]
        operation: StochasticOps,
    },

    /// Quantum physics simulations
    Quantum {
        #[command(subcommand)]
        operation: QuantumOps,
    },
    /// Display version and module information
    Info,
    /// Execute a JSON request (MCP-style interface)
    Json {
        /// JSON request string or file path (use @file.json for file, or - for stdin)
        #[arg(short, long)]
        request: String,
    },
    /// Execute JSON request from stdin (for MCP integration)
    Stdin,
    /// Start MCP server using official rmcp SDK (stdio transport)
    McpServer,
}

#[derive(Subcommand)]
enum TensorOps {
    /// Solve vacuum Einstein equations
    EinsteinSolve {
        /// Coordinate system (spherical, cartesian, etc.)
        #[arg(short, long, default_value = "spherical")]
        system: String,
    },
}

#[derive(Subcommand)]
enum CalculusOps {
    /// List available calculus operations
    List,
}

#[derive(Subcommand)]
enum FluidOps {
    /// Run cavity flow simulation
    CavityFlow {
        /// Grid size
        #[arg(short, long, default_value = "64")]
        size: usize,
        /// Reynolds number
        #[arg(short, long, default_value = "100.0")]
        reynolds: f64,
    },
}

#[derive(Subcommand)]
enum SignalOps {
    /// List available signal processing operations
    List,
}

#[derive(Subcommand)]
enum StochasticOps {
    /// List available stochastic operations
    List,
}

#[derive(Subcommand)]
enum QuantumOps {
    /// List available quantum operations
    List,
}

#[tokio::main]
async fn main() {
    let cli = Cli::parse();
    if let Some(command) = cli.command {
        run_command(command).await;
    } else {
        // Show default help when no command given
        use computational_engine::help::print_main_help;
        print_main_help();
    }
}

async fn run_command(command: Commands) {
    match command {
        Commands::Help { tool, category } => {
            use computational_engine::help_auto::print_hierarchical_help;

            if let Some(tool_name) = tool {
                // Show hierarchical help for specific tool and optional category
                print_hierarchical_help(&tool_name, category.as_deref());
            } else {
                // Show all operations organized by module
                println!("Computational Engine - All Available Operations");
                println!();
                let ops = list_all_operations();
                let mut total = 0;
                for (module, operations) in ops.iter() {
                    println!("📦 {} ({} operations):", module, operations.len());
                    total += operations.len();
                    for op in operations {
                        println!("   • {}", op);
                    }
                    println!();
                }
                println!("Total: {} operations across {} modules", total, ops.len());
                println!();
                println!("💡 TIP: For tool-specific help with examples, use:");
                println!("   brainwires-compute-engine help <TOOL> [CATEGORY]");
                println!();
                println!("Examples:");
                println!("   brainwires-compute-engine help compute           # All compute categories");
                println!("   brainwires-compute-engine help compute physics   # Physics operations only");
                println!("   brainwires-compute-engine help compute matrix    # Matrix operations only");
                println!();
                println!("Available tools: solve, differentiate, integrate, analyze, simulate,");
                println!("                 compute, transform, fieldtheory, sample, optimize");
            }
        }
        Commands::McpServer => {
            eprintln!("Starting MCP server...");
            if let Err(e) = mcp_server::server::serve_stdio().await {
                eprintln!("MCP server error: {}", e);
                std::process::exit(1);
            }
            return;
        }
        Commands::Tensor { operation } => match operation {
            TensorOps::EinsteinSolve { system } => {
                eprintln!(
                    "Solving vacuum Einstein equations in {} coordinates...",
                    system
                );
                let coords = vec![
                    "t".to_string(),
                    "r".to_string(),
                    "theta".to_string(),
                    "phi".to_string(),
                ];
                match tensor_calculus::solve_vacuum_einstein_equations(&coords, &system, &[]) {
                    Ok(solutions) => {
                        eprintln!("Found {} solutions:", solutions.len());
                        for (i, sol) in solutions.iter().enumerate() {
                            eprintln!("  Solution {}: {}", i + 1, sol.solution_type);
                        }
                    }
                    Err(e) => eprintln!("Error: {}", e),
                }
            }
        },
        Commands::Calculus { operation } => match operation {
            CalculusOps::List => {
                eprintln!("Available calculus operations:");
                eprintln!("  - Fractional derivatives and integrals");
                eprintln!("  - Special functions (Riemann zeta, elliptic integrals, etc.)");
                eprintln!("  - Variational calculus (Euler-Lagrange)");
                eprintln!("  - Stochastic calculus (Itô/Stratonovich integrals)");
                eprintln!("  - Symbolic integration");
            }
        },
        Commands::Fluid { operation } => match operation {
            FluidOps::CavityFlow { size, reynolds } => {
                eprintln!(
                    "Running cavity flow simulation (size={}, Re={})...",
                    size, reynolds
                );
                eprintln!("Simulation would run here with fluid_dynamics module");
            }
        },
        Commands::Signal { operation } => match operation {
            SignalOps::List => {
                eprintln!("Available signal processing operations:");
                eprintln!("  - FFT and inverse FFT");
                eprintln!("  - Digital filters (low-pass, high-pass, band-pass)");
                eprintln!("  - Spectrogram generation");
                eprintln!("  - Power spectral density estimation");
            }
        },
        Commands::Stochastic { operation } => match operation {
            StochasticOps::List => {
                eprintln!("Available stochastic operations:");
                eprintln!("  - Random walks and Brownian motion");
                eprintln!("  - Monte Carlo simulations");
                eprintln!("  - Markov chains");
                eprintln!("  - Stochastic differential equations");
            }
        },

        Commands::Quantum { operation } => match operation {
            QuantumOps::List => {
                eprintln!("Available quantum operations:");
                eprintln!("  - Quantum state evolution");
                eprintln!("  - Particle physics simulations");
                eprintln!("  - Tensor-based quantum computations");
                eprintln!("  - Symbolic quantum mathematics");
            }
        },
        Commands::Info => {
            println!("Computational Engine v{}", VERSION);
            println!("\nAvailable modules:");
            println!("  ✓ Tensor Calculus");
            println!("  ✓ Advanced Calculus");
            println!("  ✓ Fluid Dynamics");
            println!("  ✓ Signal Processing");
            println!("  ✓ Stochastic Processes");
            println!("  ✓ Cryptographic Mathematics");
            println!("  ✓ Symbolic Regression");
            println!("  ✓ Dimensional Analysis");
            println!("  ✓ Equation Validation");
            println!("  ✓ Computational Geometry");
            println!("  ✓ Quantum Physics (with GPU support)");
        }
        Commands::Json { request } => {
            // Handle file input if request starts with @
            let json_str = if request.starts_with('@') {
                let file_path = &request[1..];
                match std::fs::read_to_string(file_path) {
                    Ok(content) => content,
                    Err(e) => {
                        eprintln!("Error reading file '{}': {}", file_path, e);
                        std::process::exit(1);
                    }
                }
            } else {
                request.clone()
            };

            // Try new unified API first
            match serde_json::from_str::<ToolRequest>(&json_str) {
                Ok(tool_request) => {
                    // Use new unified 10-tool API
                    let dispatcher = create_default_dispatcher();
                    match dispatcher.dispatch(tool_request) {
                        Ok(response) => {
                            let output = serde_json::to_string_pretty(&response).unwrap();
                            println!("{}", output);
                        }
                        Err(e) => {
                            let error = json!({
                                "success": false,
                                "error": e
                            });
                            println!("{}", serde_json::to_string_pretty(&error).unwrap());
                            std::process::exit(1);
                        }
                    }
                }
                Err(parse_error) => {
                    let error = json!({
                        "success": false,
                        "error": format!("Invalid request format: {}. Expected ToolRequest with 'tool' and 'input' fields.", parse_error)
                    });
                    println!("{}", serde_json::to_string_pretty(&error).unwrap());
                    std::process::exit(1);
                }
            }
        }
        Commands::Stdin => {
            // Read JSON from stdin (for MCP integration)
            use std::io::{self, Read};
            let mut buffer = String::new();
            io::stdin()
                .read_to_string(&mut buffer)
                .expect("Failed to read stdin");

            // Try new unified API first
            match serde_json::from_str::<ToolRequest>(&buffer) {
                Ok(tool_request) => {
                    let dispatcher = create_default_dispatcher();
                    match dispatcher.dispatch(tool_request) {
                        Ok(response) => {
                            println!("{}", serde_json::to_string_pretty(&response).unwrap());
                        }
                        Err(e) => {
                            let error = json!({
                                "success": false,
                                "error": e
                            });
                            println!("{}", serde_json::to_string_pretty(&error).unwrap());
                            std::process::exit(1);
                        }
                    }
                }
                Err(parse_error) => {
                    let error = json!({
                        "success": false,
                        "error": format!("Invalid request format: {}. Expected ToolRequest with 'tool' and 'input' fields.", parse_error)
                    });
                    eprintln!("Failed to parse as ToolRequest: {}", parse_error);
                    println!("{}", serde_json::to_string_pretty(&error).unwrap());
                    std::process::exit(1);
                }
            }
        }
    }
}
