"""
Rename sequence IDs in .fa files based on a specified rule file.

Usage:
    python rename_seqid.py -r rule.txt -d /path/to/dir

Options:
    -r, --rule       The file containing rename rules.
    -d, --directory  The directory containing files to rename.
    -o, --output     The output directory for renamed files. Default is "rename".
    -h, --help       Show this help message and exit.

Description:
    This script renames sequence IDs in .align files based on a specified rule file.
    The rule file should contain lines in the format:
        old_id \t new_id
    Lines starting with '#' are treated as comments and ignored.

    The script will create a new directory named "rename" (or the specified output directory)
    and save the renamed files there. The renamed files will have the same names as the original files.
"""

import os
import glob
import re
import argparse
from pathlib import Path
from typing import List, Tuple
from tqdm import tqdm  # 需要安装 tqdm 库

class Color:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def validate_path(path: Path, is_dir: bool = False, check_exists: bool = True) -> None:
    """Validate input path parameters."""
    if check_exists and not path.exists():
        raise FileNotFoundError(f"{'Directory' if is_dir else 'File'} not found: {path}")
    if is_dir and not path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    if not is_dir and path.is_dir():
        raise IsADirectoryError(f"Expected file but got directory: {path}")

def load_rules(rule_file: Path) -> List[Tuple[re.Pattern, str]]:
    """Load and validate renaming rules from file."""
    compiled_rules = []
    try:
        with rule_file.open("r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                if "\t" not in line:
                    print(f"{Color.WARNING}Warning: Invalid rule format at line {line_num}: {line}{Color.ENDC}")
                    continue

                old_id, new_id = map(str.strip, line.split("\t", 1))
                if not old_id or not new_id:
                    print(f"{Color.WARNING}Warning: Empty ID at line {line_num}: {line}{Color.ENDC}")
                    continue

                try:
                    # 安全处理正则表达式特殊字符
                    pattern = re.compile(rf"^>{re.escape(old_id)}")
                    compiled_rules.append((pattern, new_id))
                except re.error as e:
                    print(f"{Color.FAIL}Error: Invalid regex pattern at line {line_num}: {e}{Color.ENDC}")
    
    except UnicodeDecodeError:
        print(f"{Color.FAIL}Error: Rule file contains invalid characters{Color.ENDC}")
        raise

    if not compiled_rules:
        raise ValueError("No valid rules found in the rule file")
    return compiled_rules

def process_alignment(input_file: Path, output_file: Path, rules: List[Tuple[re.Pattern, str]]) -> None:
    """Process a single alignment file with given rules."""
    try:
        with input_file.open("r") as fin, output_file.open("w") as fout:
            for line in fin:
                if line.startswith(">"):
                    # 应用所有规则进行匹配，将原序列ID替换为新的序列ID
                    for pattern, replacement in rules:
                        if pattern.match(line):
                            line = f">{replacement}\n"
                fout.write(line)
    except IOError as e:
        print(f"{Color.FAIL}Error processing {input_file}: {e}{Color.ENDC}")
        raise

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-r", "--rule", required=True, help="Path to the rule file")
    parser.add_argument("-d", "--directory", default=".", help="Input directory containing .fa files")
    parser.add_argument("-o", "--output", default="rename", help="Output directory for renamed files")
    
    args = parser.parse_args()

    try:
        # 验证路径参数
        rule_path = Path(args.rule)
        input_dir = Path(args.directory)
        output_dir = Path(args.output)

        validate_path(rule_path, is_dir=False)
        validate_path(input_dir, is_dir=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 加载规则
        rules = load_rules(rule_path)
        print(f"{Color.OKGREEN}Loaded {len(rules)} valid renaming rules{Color.ENDC}")

        # 获取输入文件
        input_files = list(input_dir.glob("*.align"))
        if not input_files:
            print(f"{Color.WARNING}No .fa files found in {input_dir}{Color.ENDC}")
            return

        # 处理文件
        print(f"{Color.OKBLUE}Processing {len(input_files)} files...{Color.ENDC}")
        for input_file in tqdm(input_files, desc="Processing", unit="file"):
            output_file = output_dir / input_file.name
            try:
                process_alignment(input_file, output_file, rules)
            except Exception as e:
                print(f"{Color.FAIL}Failed to process {input_file}: {e}{Color.ENDC}")

        print(f"{Color.OKGREEN}Successfully processed {len(input_files)} files{Color.ENDC}")
        print(f"Output directory: {Color.BOLD}{output_dir.resolve()}{Color.ENDC}")

    except Exception as e:
        print(f"{Color.FAIL}Error: {e}{Color.ENDC}")
        exit(1)

if __name__ == "__main__":
    main()