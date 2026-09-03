//! Shell completion generation from the live clap command graph.

use anyhow::Result;
use clap::{CommandFactory, Parser, ValueEnum};
use clap_complete::{generate, shells};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Shell {
    Bash,
    Zsh,
    Fish,
}

#[derive(Parser, Debug)]
pub struct Args {
    /// Shell whose completion script should be written to standard output.
    #[arg(value_enum)]
    pub shell: Shell,
}

/// Generate from `Cli::command()` at runtime so newly added commands cannot leave a checked-in
/// completion script stale.
pub fn run<C: CommandFactory>(args: Args) -> Result<()> {
    let mut command = C::command();
    let binary = command.get_name().to_owned();
    let mut stdout = std::io::stdout().lock();
    match args.shell {
        Shell::Bash => generate(shells::Bash, &mut command, &binary, &mut stdout),
        Shell::Zsh => generate(shells::Zsh, &mut command, &binary, &mut stdout),
        Shell::Fish => generate(shells::Fish, &mut command, &binary, &mut stdout),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{Parser, Subcommand};

    #[derive(Parser)]
    #[command(name = "fixture")]
    struct FixtureCli {
        #[command(subcommand)]
        command: FixtureCommand,
    }

    #[derive(Subcommand)]
    enum FixtureCommand {
        Run { input: String },
    }

    #[test]
    fn supported_shell_names_are_stable() {
        assert_eq!(Shell::value_variants().len(), 3);
        assert_eq!(
            Shell::from_str("bash", false).unwrap() as u8,
            Shell::Bash as u8
        );
        assert!(Shell::from_str("powershell", false).is_err());
    }

    #[test]
    fn command_factory_contains_fixture_subcommand() {
        let command = FixtureCli::command();
        assert!(command
            .get_subcommands()
            .any(|subcommand| subcommand.get_name() == "run"));
    }
}
