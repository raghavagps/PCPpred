import sys
import argparse
import pandas as pd
import re
import os
import warnings

warnings.filterwarnings('ignore')

from map_to_smiles import get_smi_from_map

# Load the MAP monomer library
df1 = pd.read_csv('data/MAP_momomers_library_new.csv')
map_to_helm_dict = df1.set_index('MAP_denotion')['Symbol'].sort_index(ascending=False).to_dict()

def smiles_is_possible(map_sequence: str) -> bool:
    """
    Returns True if a valid SMILES string can be generated, False otherwise.
    """
    smi = get_smi_from_map(map_sequence)
    # Your SMILES code returns an error message string when it fails.
    if not smi or "not possible" in smi.lower() or "cyclization not possible" in smi.lower():
        return False
    return True

##MAP to HELM sequence
def process_HELM_seq(helm_seq, ID):
    if '{' in helm_seq:
        start = helm_seq.index('{')
        end = helm_seq.index('}')
        cyc_seq = helm_seq[start:end]
        seq_len = len(helm_seq[end+1:].split('.'))
        cyc_list = cyc_seq.split('-')
        start_pos = cyc_list[0][-1]
        end_pos = cyc_list[1]
       
        if start_pos == 'N' and end_pos == 'C':
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},1:R1-{seq_len}:R2$$$'
        elif start_pos != '1' and end_pos == str(seq_len):
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},{start_pos}:R3-{seq_len}:R2$$$'
        elif start_pos == '1' and end_pos != str(seq_len):
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},1:R1-{end_pos}:R3$$$'
        else:
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},{start_pos}:R3-{end_pos}:R3$$$'
    else:
        return f'PEPTIDE{ID}{{{helm_seq}}}$$$$'
def convert_map_to_helm_sequence(map_str, ID):
    nterm_pattern = r'\{nt:[^}]+\}'
    cyc_pattern = r'\{cyc:\s*([N]|\d+)-([C]|\d+)\}'
    string = ''
    nterm_modifications = re.findall(nterm_pattern, map_str)
    map_str = re.sub(nterm_pattern, '', map_str)
    cyc_string = re.search(cyc_pattern, map_str)
    # print("cyc_string",cyc_string[0])
    map_str = re.sub(cyc_pattern, '', map_str)
    if cyc_string:
        string += ''.join(cyc_string[0]) + ''.join(nterm_modifications) + map_str
    else:
        string += ''.join(nterm_modifications) + map_str

    tokens = []
    i = 0
    while i < len(string):
        matched = False
        for key in map_to_helm_dict.keys():
            if string[i:].startswith(key):
                if string[i-4:i] == 'cyc:':
                    val = map_to_helm_dict[key]
                    token = f'[{val}]' if len(val) > 1 else f'{val}'
                    tokens.append(token)
                    i += len(key)
                    matched = True
                    break
                else:
                    val = map_to_helm_dict[key]
                    token = f'[{val}].' if len(val) > 1 else f'{val}.'
                    tokens.append(token)
                    i += len(key)
                    matched = True
                    break

        if not matched:
            tokens.append(string[i])
            i += 1
    helm_seq = ''.join(tokens).rstrip('.')
    # print(helm_seq)
    return helm_seq

import datetime
def main():
    parser = argparse.ArgumentParser(description='Convert MAP sequence(s) to HELM sequence(s)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--sequence', type=str,
                       help='Single FASTA-like input: header line starting with ">" followed by MAP sequence line (use "\\n" for newline)')
    group.add_argument('-f', '--file', type=str,
                       help='Input file with FASTA-like MAP sequences')
    parser.add_argument('-id', '--identifier', type=str,
                       help='ID used in HELM notation (required with -s, ignored with -f)')
    parser.add_argument('-o', '--output', type=str,
                       help='Output file path for HELM (used with -f)')

    args = parser.parse_args()

    if args.sequence:
        if not args.identifier:
            parser.error("--id/--identifier is required when using --sequence")
        ID = args.identifier

        # Handle escaped \n in input
        seq_input = args.sequence.encode().decode('unicode_escape')
        parts = seq_input.strip().split("\n")
        if len(parts) < 2 or not parts[0].startswith(">"):
            print("Error: Input must be FASTA-like with a header starting with '>' and a MAP sequence line.")
            sys.exit(1)
        header = parts[0][1:].strip()
        seq = parts[1].strip()

        if smiles_is_possible(seq):
            helm_seq = convert_map_to_helm_sequence(seq, ID)
            result = process_HELM_seq(helm_seq, ID)
            print(f">{header}\n{result}")
        else:
            print(f">{header}\nHELM conversion not possible (no valid SMILES)")

    elif args.file:
       
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = args.output if args.output else os.path.join(
            'results', f"helm_output_{timestamp}.txt"
        )
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)

        with open(args.file, 'r') as f, open(output_file, 'w') as out:
            header, seq_line = None, None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header and seq_line:
                        try:
                            map_seq, peptide_id = seq_line.split(',')
                            if smiles_is_possible(map_seq.strip()):
                                helm_seq = convert_map_to_helm_sequence(map_seq.strip(), peptide_id.strip())
                                result = process_HELM_seq(helm_seq, peptide_id.strip())
                                out.write(f"{header}\t{map_seq.strip()}\t{result}\n")
                            else:
                                out.write(f"{header}\t{map_seq.strip()}\tHELM conversion not possible (no valid SMILES)\n")
                        except ValueError:
                            out.write(f"Error: Invalid format in line '{seq_line}'. Expected 'MAP_sequence,ID'\n")
                    header = line[1:].strip()
                    seq_line = None
                else:
                    seq_line = line
            # Process the last record
            if header and seq_line:
                try:
                    map_seq, peptide_id = seq_line.split(',')
                    if smiles_is_possible(map_seq.strip()):
                        helm_seq = convert_map_to_helm_sequence(map_seq.strip(), peptide_id.strip())
                        result = process_HELM_seq(helm_seq, peptide_id.strip())
                        out.write(f"{header}\t{map_seq.strip()}\t{result}\n")
                    else:
                        out.write(f"{header}\t{map_seq.strip()}\tHELM conversion not possible (no valid SMILES)\n")
                except ValueError:
                    out.write(f"Error: Invalid format in line '{seq_line}'. Expected 'MAP_sequence,ID'\n")

        print(f"HELM sequences written to {output_file}")


if __name__ == "__main__":
    main()
