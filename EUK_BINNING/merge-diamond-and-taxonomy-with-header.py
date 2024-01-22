#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

def main(args):
    if not args.diamond_output:
        raise ConfigError("You must provide a DIAMOND output file path.")

    if not args.taxonomy:
        raise ConfigError("You must provide a taxonomy file path.")

    filesnpaths.is_file_tab_delimited(args.diamond_output)
    filesnpaths.is_file_tab_delimited(args.taxonomy)

    progress.new("Bleep bloop")
    progress.update('...')
    diamond = [l.strip().split('\t') for l in open(args.diamond_output).readlines()]
    taxonomy = utils.get_TAB_delimited_file_as_dictionary(args.taxonomy)
    fields_in_taxonomy_file = utils.get_columns_of_TAB_delim_file(args.taxonomy)
    
    with open(args.out, 'w') as output:
        output.write("Query"+'\t'+"Target"+'\t'+"Percent Identity"+'\t'+"Alignment length"+'\t'+"Mismatches"+'\t'+"gaps"+'\t'+"start aligment query"+'\t'+"end alignment query"+'\t'+"start aligment target"+'\t'+"end alignment target"+'\t'+"evalue"+'\t'+"bitscore"+'\t'+"Paralog_Type"+'\t'+"Source"+'\t'+"Domain"+'\t'+"Super_Group"+'\t'+"Phylum"+'\t'+"Class"+'\t'+"Order"+'\t'+"Family"+'\t'+"Genre"+'\n')
        for line in diamond:
            if line[1] in taxonomy:
                for field in fields_in_taxonomy_file[1:]:
                    if taxonomy[line[1]][field]:
                        line.append(taxonomy[line[1]][field])
                    else:
                        line.append('')
            else:
                for field in fields_in_taxonomy_file[1:]:
                    line.append('')

            output.write('\t'.join(line) + '\n')
    progress.end()

    run.info('Merged output', args.out)


if __name__ == '__main__':
    parser = ArgumentParser(description="Expand diamond output with taxonomy. The first column of the taxonomy file must match to "
                "the second column of the DIAMOND output file.")
    parser.add_argument('--diamond-output', metavar = 'FILE',
                                        help = 'DIAMOND search output')
    parser.add_argument('--taxonomy', metavar = 'FILE',
                                        help = 'Taxonomy file')
    parser.add_argument('--out', metavar = 'FILE',
                                        help = 'output file')

    args = parser.parse_args()

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
