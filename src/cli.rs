//! CLI argument parsing for configuring what to do

use crate::{AppSettings, ChildDirection, Float, TreePath};
use std::ffi::OsString;
use std::fmt::Display;
use std::path::PathBuf;
use std::process::exit;
use std::str::FromStr;

/// How should we give information about the model back to the user?
pub enum DisplayMethod<'a> {
    /// Output data about the model in a CSV format, printing directly to Stdout
    Csv {
        file: Option<&'a str>,
        add_paths: Vec<TreePath>,
    },
    /// Produce a series of PNG files, each named 'branch-tree-*.png' in the current directory.
    Png { file_pattern: &'a str },
}

/// Which model of the lungs should we use?
pub enum Model {
    /// `EqualChildGenerator` with the given depth
    Symmetric { depth: usize },
    /// `FromJsonGenerator`, loading from the given file
    FromJson { file: PathBuf },
}

/// Main entrypoint for the program
///
/// Handles command-line argument parsing & others, then dispatches based on what's been enabled.
pub fn run() {
    let matches = args_parser().get_matches();

    let total_time = matches
        .value_of("total-time")
        .unwrap()
        .parse::<Float>()
        .expect("validator should ensure that this is a valid float");

    let timestep = matches
        .value_of("timestep")
        .unwrap()
        .parse::<Float>()
        .expect("validator should ensure that this is a valid float");

    let display_method = match matches.value_of("output").unwrap() {
        "csv" | "CSV" => {
            let mut add_paths = Vec::new();
            if let Some(values) = matches.values_of("add-csv-info") {
                for p in values {
                    match parse_tree_path(p) {
                        Ok(p) => add_paths.push(p),
                        Err(e) => {
                            eprintln!("invalid tree path {:?}: {}", p, e);
                            exit(1)
                        }
                    }
                }
            }

            DisplayMethod::Csv {
                file: matches.value_of("file-pattern"),
                add_paths,
            }
        }
        "png" | "PNG" => {
            // 'file-pattern' is required for 'png' in the clap app
            let pat = matches.value_of("file-pattern").unwrap();

            if matches.is_present("add-csv-info") {
                eprintln!("-a/--add-csv is invalid in PNG output mode");
                exit(1);
            }

            if pat.matches("{}").count() != 1 {
                eprintln!(
                    "PNG file pattern must contain exactly one occurence of '{{}}' to substitute"
                );
                exit(1);
            } else {
                DisplayMethod::Png { file_pattern: pat }
            }
        }
        _ => unreachable!(),
    };

    let model = match matches.value_of("lung-model").unwrap() {
        "symmetric" => {
            let depth = matches
                .value_of("depth")
                .unwrap()
                .parse::<usize>()
                .expect("validator should ensure this is a valid usize");

            Model::Symmetric { depth }
        }
        "from-json" => {
            let file: OsString = matches.value_of_os("input-file").unwrap().to_owned();

            Model::FromJson {
                file: PathBuf::from(file),
            }
        }
        _ => unreachable!(),
    };

    AppSettings {
        total_time,
        timestep,
        display_method,
        model,
    }
    .run()
}

fn parse_tree_path(mut s: &str) -> Result<TreePath, String> {
    let mut path = Vec::new();

    while !s.is_empty() {
        if s.starts_with(".left") {
            path.push(ChildDirection::Left);
            s = &s[".left".len()..];
        } else if s.starts_with(".right") {
            path.push(ChildDirection::Right);
            s = &s[".right".len()..];
        } else {
            return Err(format!(r#"expected ".left" or ".right", found {:?}"#, s));
        }
    }

    if path.is_empty() {
        return Err("path must be non-empty".to_owned());
    }

    Ok(path)
}

/// Returns the String generated by substituting the image number into a PNG file pattern
///
/// Panics if the pattern is invalid. Any pattern given to an `AppSettings` by this module will be
/// valid.
///
/// The image number also does not *necessarily* have to be less than `num_imgs` -- that's actually
/// the case when we're given total_time = 0, to repeat forever.
pub fn substitute_png_file_pattern(pat: &str, img_number: u64, mut num_imgs: u64) -> String {
    let mut digits_space = 0;
    while num_imgs > 0 {
        num_imgs /= 10;
        digits_space += 1;
    }

    let (fst_half, snd_half) = pat.split_once("{}").expect("pattern should contain '{}'");
    format!(
        "{}{:0digits$}{}",
        //  ^^^^^^^^^ pad with leading zeroes so there's at least 'digits' digits
        fst_half,
        img_number,
        snd_half,
        digits = digits_space
    )
}

fn args_parser() -> clap::App<'static, 'static> {
    use clap::{App, Arg};

    // Generic validator, to check that we can parse the value as a certain type.
    fn valid_as<T>(check: fn(T) -> bool, msg: &'static str) -> impl Fn(String) -> Result<(), String>
    where
        T: FromStr,
        T::Err: Display,
    {
        move |s| {
            s.parse::<T>()
                .map_err(|e| e.to_string())
                .and_then(|v| match check(v) {
                    false => Err(msg.to_owned()),
                    true => Ok(()),
                })
        }
    }

    App::new("Super-fancy Airflow Simulator")
        .version("0.1.0")
        .author("Max Sharnoff <github@max.sharnoff.org>")
        .about("Tool for simulating airflow in the lungs")
        ////////////////////////////////////////////////////////////
        // SECTION: Simulation meta-variables                     //
        ////////////////////////////////////////////////////////////
        .arg(
            Arg::with_name("total-time")
                .short("t")
                .long("time")
                .value_name("SECONDS")
                .help("Set how long to run for, in seconds. Use 0 to run until stopped")
                .takes_value(true)
                .required(true)
                .validator(valid_as::<Float>(|f| f >= 0.0, "total time must be >= 0")),
        )
        .arg(
            Arg::with_name("timestep")
                .long("timestep")
                .value_name("SECONDS")
                .help("Set the timestep to advance the simulation by each time, in seconds")
                .long_help(concat!(
                    "Set the timestep to advance the simulation by each time, in seconds.\n",
                    "Defaults to 0.01. Also determines the frequency of PNG output images.",
                ))
                .default_value("0.01")
                .takes_value(true)
                .required(true)
                .validator(valid_as::<Float>(|f| f > 0.0, "timestep must be > 0")),
        )
        ////////////////////////////////////////////////////////////
        // SECTION: Output methods                                //
        ////////////////////////////////////////////////////////////
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .possible_values(&["csv", "CSV", "png", "PNG"])
                .requires_ifs(&[
                    // 'file-pattern' is required for emitting PNGs
                    ("png", "file-pattern"),
                    ("PNG", "file-pattern"),
                ])
                .required(true)
                .takes_value(true)
                .help("Set the output format to use: 'csv' or 'png'")
                .long_help(concat!(
                    "Set the output format to use: 'csv' or 'png'. The CSV output\n",
                    "is an overview of the airflow & total volume at each timestep.\n",
                    "PNG output provides an image at each timestep, and requires the\n",
                    "'-f'/'--file-pattern' option.\n\n",
                    "CSV output defaults to Stdout if no file is provided",
                )),
        )
        .arg(
            Arg::with_name("file-pattern")
                .short("f")
                .long("file-pattern")
                .value_name("FILE-PATTERN")
                .help("Set the file pattern to use when outputting PNGs or a single CSV")
                .long_help(concat!(
                    "Set the file pattern to use when outputting PNGs or a single CSV.\n",
                    "When making PNGs, the pattern should have a literal '{}' somewhere\n",
                    "which will be replaced with an incrementing number for the image.",
                )),
        )
        .arg(
            Arg::with_name("add-csv-info")
                .short("a")
                .long("add-csv")
                .value_name("TREE-PATH")
                .multiple(true)
                .help("Add a path in the tree to include in the CSV output")
                .long_help(concat!(
                    "Add a path in the tree to include in the CSV output.\n",
                    "For a given $path, this will add the fields '$path flow out' and\n",
                    "'$path total volume'.\n",
                    "The path must be written as .(left|right).(left|right)..., for example:\n",
                    "\n",
                    "    .left.right.left    or    .right\n",
                    "\n",
                    "are both valid paths.",
                )),
        )
        ////////////////////////////////////////////////////////////
        // SECTION: Model selection & customization               //
        ////////////////////////////////////////////////////////////
        .arg(
            Arg::with_name("lung-model")
                .short("m")
                .long("model")
                .possible_values(&["symmetric", "from-json"])
                .requires_ifs(&[("symmetric", "depth"), ("from-json", "input-file")])
                .required(true)
                .takes_value(true)
                .help("Set the model of the lungs to use: 'symmetric' or 'from-json'")
                .long_help(concat!(
                    "Set the model of the lungs to use: current options are 'symmetric' or\n",
                    "'from-json'\n\n",
                    "'symmetric' requires the '--depth' parameter to be given, which determines\n",
                    "what the terminal child depth should be.\n\n",
                    "'from-json' requires the '--input-file' parameter to be given",
                )),
        )
        .arg(
            Arg::with_name("depth")
                .long("depth")
                .value_name("DEPTH")
                .help("Set the terminal child depth of the symmetric model")
                .takes_value(true)
                .validator(valid_as::<usize>(|d| d > 0, "depth must be > 0")),
        )
        .arg(
            Arg::with_name("input-file")
                .long("input-file")
                .value_name("FILE")
                .help("Sets the input file to load a JSON model from")
                .takes_value(true),
        )
}
