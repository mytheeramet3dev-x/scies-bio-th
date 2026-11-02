#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioError {
    InvalidCharacter(char),
    ParseError(String),
    AlignmentMismatch,
    TranslationError,
    ReprogrammingFailed(String),
    IncompatibleGamete(String),
    DevelopmentArrested(String),
    SimulationError(String),
    InvalidInput(String),
    IoError(String),
}

impl std::fmt::Display for BioError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BioError::InvalidCharacter(c) => write!(f, "Invalid character: {}", c),
            BioError::ParseError(e) => write!(f, "Parse error: {}", e),
            BioError::AlignmentMismatch => write!(f, "Alignment mismatch"),
            BioError::TranslationError => write!(f, "Translation error"),
            BioError::ReprogrammingFailed(e) => write!(f, "Epigenetic reprogramming failed: {}", e),
            BioError::IncompatibleGamete(e) => write!(f, "Incompatible gametes: {}", e),
            BioError::DevelopmentArrested(e) => write!(f, "Embryonic development arrested: {}", e),
            BioError::SimulationError(e) => write!(f, "Simulation error: {}", e),
            BioError::InvalidInput(e) => write!(f, "Invalid input: {}", e),
            BioError::IoError(e) => write!(f, "I/O error: {}", e),
        }
    }
}

impl std::error::Error for BioError {}

pub type Result<T> = std::result::Result<T, BioError>;
