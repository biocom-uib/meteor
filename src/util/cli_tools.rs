use std::{process::{Command, Stdio}, path::{Path, PathBuf}, io::{ErrorKind, self}};

pub trait CliTool {
    fn new_command(&self) -> Command;
}

#[derive(Debug, Clone)]
pub struct BlastTool {
    install_prefix: Option<PathBuf>,
    command: String,
}

pub fn silent_test_command_exists(mut command: Command) -> io::Result<bool> {
    let result = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();

    match result {
        Err(e) => {
            if e.kind() == ErrorKind::NotFound {
                Ok(false)
            } else {
                Err(e)
            }
        },
        Ok(_) => {
            Ok(true)
        }
    }
}

pub fn resolve_with_silent_test<Tool, CmdBuilder>(
    tool: Tool,
    builder: CmdBuilder,
) -> io::Result<Option<Tool>>
where
    Tool: CliTool,
    CmdBuilder: FnOnce(&mut Command),
{
    let mut test_command = tool.new_command();

    builder(&mut test_command);

    if silent_test_command_exists(test_command)? {
        Ok(Some(tool))
    } else {
        Ok(None)
    }
}

pub fn resolve_with_silent_test_args<Tool>(tool: Tool, args: &[&str]) -> io::Result<Option<Tool>>
where
    Tool: CliTool,
{
    resolve_with_silent_test(tool, |cmd| {
        cmd.args(args);
    })
}

impl BlastTool {
    pub fn resolve(install_prefix: Option<&Path>, command: &str) -> io::Result<Option<BlastTool>> {
        let blast_tool = BlastTool {
            install_prefix: install_prefix.map(Path::to_owned),
            command: command.to_owned(),
        };

        resolve_with_silent_test_args(blast_tool, &["-version"])
    }
}

impl CliTool for BlastTool {
    fn new_command(&self) -> Command {
        if let Some(prefix) = &self.install_prefix {
            Command::new(prefix.join("bin").join(&self.command))
        } else {
            Command::new(&self.command)
        }
    }
}

#[derive(Debug, Clone)]
pub struct Java;

impl Java {
    pub fn resolve() -> io::Result<Option<Self>> {
        resolve_with_silent_test_args(Java {}, &["-version"])
    }
}

impl CliTool for Java {
    fn new_command(&self) -> Command {
        Command::new("java")
    }
}


#[derive(Debug, Clone)]
pub enum Prodigal {
    InPath(String),
    SpecifiedPath(PathBuf),
}

impl Prodigal {
    pub fn resolve_at(path: &Path) -> io::Result<Option<Self>> {
        let prodigal_tool = Prodigal::SpecifiedPath(path.to_path_buf());

        resolve_with_silent_test_args(prodigal_tool, &["-v"])
    }

    fn resolve_name_in_path(name: &str) -> io::Result<Option<Self>> {
        let prodigal_tool = Prodigal::InPath(name.to_string());

        resolve_with_silent_test_args(prodigal_tool, &["-v"])
    }

    pub fn resolve_prodigal_in_path() -> io::Result<Option<Self>> {
        if let result@Some(_) = Self::resolve_name_in_path("prodigal")? {
            Ok(result)
        } else if let result@Some(_) = Self::resolve_name_in_path("prodigal.linux")? {
            Ok(result)
        } else {
            Ok(None)
        }
    }

    pub fn resolve(path: Option<&Path>) -> io::Result<Option<Self>> {
        if let Some(path) = path {
            Self::resolve_at(path)
        } else {
            Self::resolve_prodigal_in_path()
        }
    }
}

impl CliTool for Prodigal {
    fn new_command(&self) -> Command {
        match self {
            Prodigal::InPath(command) => Command::new(command),
            Prodigal::SpecifiedPath(path) => Command::new(path),
        }
    }
}
