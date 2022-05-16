#!/bin/python
# 
# Helper script for filling filling out pgfplots with data derived from the 'data' directory. It's
# heavily inspired by a script of mine that converts CSVs to a series of tuples compatable with
# pgfplots -- so the choices here reflect that.

import argparse
import csv
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", nargs=1,
            help="Sets the input file to use. Must end in .autoplot.tex")

    args = parser.parse_args()
    infile: str = args.inputfile[0]
    if not infile.endswith('.autoplot.tex'):
        raise Exception(f"Input file must end in '.autoplot.tex'. Got: {infile}")

    # outfile = inputfile : s/.autoplot.tex$/.generated.tex/
    outfile = infile[:len(infile) - len('.autoplot.tex')] + '.generated.tex'
    
    # Run:
    run_for_file(infile, outfile)

def run_for_file(infile: str, outfile: str) -> None:
    marker = r'% ?@generated\b'
    pattern = re.compile(r'^(\s*)'+marker+r'(.*)$')

    setfileprefix = re.compile(r'^\s*% ?@setfileprefix(\((?P<prefix>.*)\)$)?')

    try:
        inf = open(infile, 'r')
        outf = open(outfile, 'w+')

        context = {}

        for line_idx, line in enumerate(inf.readlines()):
            stripped = line.rstrip('\r\n')
            line_end = line[len(stripped):]

            if (m := pattern.fullmatch(stripped)) is not None:
                leading_whitespace = m.group(1)
                lineno = line_idx + 1
                args = GenerateArgs(m.group(2), lineno) # line number provided for context
                line = leading_whitespace + process_line(context, args, lineno) + line_end
            elif (m := setfileprefix.fullmatch(stripped)) is not None:
                if (p := m.group('prefix')) is None:
                    raise Exception(f'line {line_idx + 1}: malformed @setfileprefix')
                context['prefix'] = p
                continue # skip this line when writing output
            
            outf.write(line)

    finally:
        if 'inf' in locals():
            locals()['inf'].close() # go through 'locals' so pyright doesn't get mad at us
        if 'outf' in locals():
            locals()['outf'].close()

class GenerateArgs:
    def __init__(self, string_description: str, lineno: int) -> None:
        paren_value = r"'(?P<quotevalue>[^']*)'"
        normal_value = r"[^,\s\(\)]+"
        value = f'(?P<value>{paren_value}|{normal_value})'

        argname = r'(?P<argname>[a-z][a-z0-9_\-]*)'
        arg = f'({argname}=)?{value}'
        
        trailing_comma = r'(?P<trailing_comma>,\s*)'
        pattern = re.compile(f'{arg}({trailing_comma}|$)')

        first_pass_parens = string_description.startswith('(') and string_description.endswith(')')
        if not first_pass_parens:
            raise Exception(f"error: line {lineno} expected @generated(...), with nothing after closing parenthesis")

        middle = string_description[1:-1]
        args = {}
        # arg name, with aliases
        required_order = [('fields', ['f', 'fs']), ('range', ['r'])]
        optional = [('stride', ['s']), ('transform', ['t', 'tx'])]

        # args are parsed in the order of: required..., maybe optional, (unnamed) filepath
        while middle != '':
            arg_match = pattern.match(middle)
            if arg_match is None:
                raise Exception(f'line {lineno}: malformed arguments to @generated')

            match_len = len(arg_match.group(0))
            name = arg_match.group('argname')
            quoteval = arg_match.group('quotevalue')
            val = arg_match.group('value') if quoteval is None else quoteval # maybe unpack quotes
            trailing = arg_match.group('trailing_comma')

            if (name is None) != (trailing is None):
                raise Exception(f'line {lineno}: either argname must be present or filename at end of arguments')

            if len(required_order) != 0:
                next_name, aliases = required_order[0]

                if name is None: # don't jump ahead to the filepath
                    raise Exception(f'line {lineno}: missing required arguments')
                elif name != next_name and name not in aliases:
                    raise Exception(f'line {lineno}: expected argname {next_name}, found {name}')

                args[next_name] = val
                del required_order[0]
            elif name is not None: # name present, must be optional arg
                for i, (n, aliases) in enumerate(optional):
                    if name == n or name in aliases:
                        args[n] = val
                        del optional[i]
                        break
                else:
                    raise Exception(f'line {lineno}: unrecognized or duplicate arg {name}')
            else: # no name, must be the final path. Regex says this will be the end of the string
                args['path'] = val

            middle = middle[match_len:]

        # done grabbing the arguments. Now gotta parse & validate
        if len(required_order) != 0:
            raise Exception(f'line {lineno}: missing required arguments')
        elif 'path' not in args:
            raise Exception(f'line {lineno}: missing final path argument')

        fieldspec = args['fields'].split(',')
        if len(fieldspec) != 2:
            raise Exception(f'line {lineno}: fields must be specified as: <fieldnum>,<fieldnum>')
        try:
            self.x_field = int(fieldspec[0])
            self.y_field = int(fieldspec[1])
        except BaseException:
            raise Exception(f'line {lineno}: fields must be integers')
        if self.x_field <= 0 or self.y_field <= 0:
            raise Exception(f'line {lineno}: fields use 1-based indexing, must be greater than zero')

        self.line_range = Range.parse(lineno, args['range'])

        if (s := args.get('stride')) is not None:
            try:
                self.stride = int(s)
            except BaseException:
                raise Exception(f'line {lineno}: stride must be an integer')
        else:
            self.stride = 1

        if (t := args.get('transform')) is not None:
            try:
                self.transform = eval(t, {}, {})
            except BaseException:
                raise Exception(f'line {lineno}: failed to evaluate transform')
        else:
            self.transform = lambda x, y: (x, y)

        self.path = args['path']

# This class taken in full from my personal script
class Range:
    def __init__(self, start: int | None, end: int | None):
        self.start = start
        self.end = end

    @staticmethod
    def parse(lineno: int, s: str) -> 'Range':
        def raise_parse_error():
            raise Exception(f'line {lineno}: range must match the format [ start ] ".." [ [ "=" ] end ]' )
        
        parts = s.split('..')
        if len(parts) != 2:
            raise_parse_error()

        start_str, end_str = parts
        if start_str == '':
            start = None
        elif start_str.isnumeric():
            start = int(start_str)
        else:
            raise_parse_error()

        if end_str == '':
            end = None
        else:
            if end_str.startswith('='):
                end_str = end_str[1:]
                end_exclusive = True
            else:
                end_exclusive = False

            if not end_str.isnumeric():
                raise_parse_error()

            end = int(end_str)
            if end_exclusive:
                end += 1
        
        return Range(start, end)

    def contains(self, lineno: int) -> bool:
        return (self.start is None or lineno >= self.start) and (self.end is None or lineno < self.end)

def process_line(context: dict[str, str], args: GenerateArgs, lineno: int) -> str:
    prefix = context.get('prefix')
    filepath = args.path if prefix is None else prefix + args.path

    # make filepath relative to git repo root, not 'paper':
    filepath = f'../{filepath}'

    # go from 1-based indexing to 0-based
    args.x_field -= 1
    args.y_field -= 1

    output = ''
    rowno = 0
    stride_counter = 1 # initially set to 1 so that we start at the first valid line

    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            rowno += 1 # pre-add to start line numbers at 1
            if not args.line_range.contains(rowno):
                if output != '': # previously in range, now out of range
                    break
                continue # not in range yet

            stride_counter -= 1
            if stride_counter != 0:
                continue
            else:
                # stride_counter = 0; reset it
                stride_counter = args.stride

            if len(row) < max(args.x_field, args.y_field):
                raise Exception(f'line {lineno}: input file does not have enough fields on line {rowno}')

            x, y = (row[args.x_field], row[args.y_field])
            new_pair = args.transform(x, y)
            if not isinstance(new_pair, (list, tuple)) or len(new_pair) != 2:
                raise Exception(f'line {lineno}: transform did not return a list or tuple with length 2')

            x, y = new_pair
            output += make_point(x, y)

    if args.line_range.start is not None and rowno < args.line_range.start:
        raise Exception(f'line {lineno}: input is not long enough: range start is {args.line_range.start} but only {rowno} lines of input')
    elif args.line_range.end is not None and rowno < (end := args.line_range.end - 1):
        raise Exception(f'line {lineno}: input is not long enough: range end is {end} but only {rowno} lines of input')

    return output

def make_point(x: str, y: str) -> str:
    return f'({x},{y})'

if __name__ == '__main__':
    main()
