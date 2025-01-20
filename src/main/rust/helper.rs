
use anyhow;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum PublicError
{
    #[error("ApplicationError: {0}")]
    ApplicationError (String),

//    #[error("DataError: {0}")]
//    DataError (String)
}

impl From<anyhow::Error> for PublicError
{
    fn from (err: anyhow::Error)
        -> PublicError
    {
        PublicError::ApplicationError (err.to_string ())
    }
}

