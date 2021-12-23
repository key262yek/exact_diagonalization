// Moduel for defining errors
use std::fmt::{self, Display, Formatter};
use std::io;


#[derive(Debug)]
pub struct Error{
    /// This `Box` allows us to keep the size of `Error` as small as possible. A
    /// larger `Error` type was substantially slower due to all the functions
    /// that pass around `Result<T, Error>`.
    err: Box<ErrorCode>,
}

pub const MAX_TRIAL : usize = 100;

/// Alias for a `Result` with the error type `serde_json::Error`.
pub type Result<T> = std::result::Result<T, Error>;

impl Error{
    pub fn make_error_msg(msg: String) -> Self{                         // message만을 담고 있는 error
        Error {
            err: Box::new(ErrorCode::Message(msg.into_boxed_str())),
        }
    }

    pub fn make_error_io(error: io::Error) -> Self{                     // io에서 돌아온 error
        Error {
            err: Box::new(ErrorCode::Io(error)),
        }
    }

    pub fn make_error_syntax(code : ErrorCode) -> Self{                 // 직접 정의한 error들. syntax error 중심
        Error {
            err: Box::new(code),
        }
    }

    pub fn classify(&self) -> Category{                                 // error들을 분류
        match *self.err{
            ErrorCode::Message(_) => Category::Data,
            ErrorCode::Io(_) => Category::Io,
            ErrorCode::InvalidBitIndex
            | ErrorCode::OverFlow
             => Category::Syntax,
        }
    }

    pub fn is_io(&self) -> bool{
        self.classify() == Category::Io
    }

    pub fn is_data(&self) -> bool{
        self.classify() == Category::Data
    }

    pub fn is_syntax(&self) -> bool{
        self.classify() == Category::Syntax
    }
}

impl PartialEq for Error{
    // 여러 이유로 Error 구조체는 등식을 정의하기 어렵다.
    // 그래서 대신 category만 같으면 같은 에러로 취급하는 것
    fn eq(&self, other : &Self) -> bool{
        self.classify() == other.classify()
    }
}


/// Categorizes the cause of a `serde_json::Error`.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Category{
    Io,
    Data,
    Syntax,
}

#[derive(Debug)]
pub enum ErrorCode{
    /// Catchall for syntax error message
    Message(Box<str>),

    /// Some IO error occurred while serializing or deserializing.
    Io(io::Error),

    // Invalid Bit index.
    InvalidBitIndex,

    // Over Flow
    OverFlow,
}

impl Display for ErrorCode{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result{
        match self{
            ErrorCode::Message(ref msg) => f.write_str(msg),
            ErrorCode::Io(ref err) => Display::fmt(err, f),
            ErrorCode::InvalidBitIndex => f.write_str("Invalid Bit index."),
            ErrorCode::OverFlow => f.write_str("State representation over flows."),
        }
    }
}

impl Display for Error{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result{
        Display::fmt(&*self.err, f)
    }
}


#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_fmt(){
        use std::io;

        assert_eq!(format!("{}", Error::make_error_msg("Test message for Message error".to_string())).as_str(),
            "Test message for Message error");
        assert_eq!(format!("{}", Error::make_error_io(io::Error::new(io::ErrorKind::NotFound,  "Test io error"))).as_str(),
            "Test io error");
        assert_eq!(format!("{}", Error::make_error_syntax(ErrorCode::InvalidBitIndex)).as_str(),
            "Invalid Bit index.");
        assert_eq!(format!("{}", Error::make_error_syntax(ErrorCode::OverFlow)).as_str(),
            "State representation over flows.");
    }

    #[test]
    fn test_classify(){
        use std::io;

        assert_eq!(Error::make_error_msg("Test message".to_string()).classify(),
            Category::Data);
        assert_eq!(Error::make_error_io(io::Error::new(io::ErrorKind::NotFound,  "Test io error")).classify(),
            Category::Io);
        assert_eq!(Error::make_error_syntax(ErrorCode::InvalidBitIndex).classify(), Category::Syntax);
        assert_eq!(Error::make_error_syntax(ErrorCode::OverFlow).classify(), Category::Syntax);
    }
}
