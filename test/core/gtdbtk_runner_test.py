import json
import logging
import os
import tempfile

from pathlib import Path

from kb_gtdbtk.core.gtdbtk_runner import run_gtdbtk

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s', level=logging.INFO)


def test_gtdbtk_run():

    with tempfile.TemporaryDirectory(prefix='test_gtdbtk_run') as test_dir_str:
        test_dir = Path(test_dir_str)
        out_dir = test_dir / 'output'
        out_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = test_dir / 'temp'
        temp_dir.mkdir(parents=True, exist_ok=True)

        tf = []

        def runner(command):
            tf.append(command.pop(5))

            td = command[3]
            assert Path(td).is_dir()

            temp_out = temp_dir / 'output'

            assert command == [
                'gtdbtk',
                'classify_wf',
                '--out_dir', str(temp_out),
                '--batchfile',  # arg popped
                '--cpus', '16',
                '--min_perc_aa', '50.2']

            # arbitrary TSV files, these may not match what GTDB-tk produces
            # implmentation doesn't care for now other than the first column

            with open(temp_out / 'gtdbtk.ar122.summary.tsv', 'w') as t:
                t.writelines(['\t'.join(['Name', 'field1', 'field2']) + '\n',
                              '\t'.join(['a1', 'foo', 'bar']) + '\n',
                              '\t'.join(['a2', 'baz', 'bat']) + '\n',
                              ])

            with open(temp_out / 'gtdbtk.bac120.summary.tsv', 'w') as t:
                t.writelines(['\t'.join(['Name', 'field1', 'field2']) + '\n',
                              '\t'.join(['a1', 'whoo', 'whee']) + '\n',
                              '\t'.join(['a2', 'whoa', 'whump']) + '\n',
                              ])

            with open(temp_out / 'gtdbtk.bac120.markers_summary.tsv', 'w') as t:
                t.writelines(['\t'.join(['user_genome', 'field1', 'field2']) + '\n',
                              '\t'.join(['a1', 'fee', 'fie']) + '\n',
                              '\t'.join(['a2', 'fo', 'fum']) + '\n',
                              ])

            # skip 'gtdbtk.ar122.markers_summary.tsv'

            # ####   end of runner callable   ####

        run_gtdbtk(
            runner,
            {'a1': {'path': '/somepath1', 'assembly_name': 'somefile1.fasta'},
             'a2': {'path': '/somepath2', 'assembly_name': 'somefile2.fasta'},
             },
            out_dir,
            temp_dir,
            50.2,
            16
            )

        with open(tf[0]) as bf:
            lines = bf.readlines()
            assert len(lines) == 2
            assert lines[0] == '/somepath1\ta1\n'
            assert lines[1] == '/somepath2\ta2\n'

        assert sorted(os.listdir(out_dir)) == [
            'gtdbtk.ar122.summary.tsv',
            'gtdbtk.ar122.summary.tsv.json',
            'gtdbtk.bac120.markers_summary.tsv',
            'gtdbtk.bac120.markers_summary.tsv.json',
            'gtdbtk.bac120.summary.tsv',
            'gtdbtk.bac120.summary.tsv.json',
        ]

        with open(out_dir / 'gtdbtk.ar122.summary.tsv.json') as j:
            assert json.load(j) == {'data': [
                {'Name': 'somefile1.fasta', 'field1': 'foo', 'field2': 'bar'},
                {'Name': 'somefile2.fasta', 'field1': 'baz', 'field2': 'bat'},
            ]}

        with open(out_dir / 'gtdbtk.bac120.summary.tsv.json') as j:
            assert json.load(j) == {'data': [
                {'Name': 'somefile1.fasta', 'field1': 'whoo', 'field2': 'whee'},
                {'Name': 'somefile2.fasta', 'field1': 'whoa', 'field2': 'whump'},
            ]}

        with open(out_dir / 'gtdbtk.bac120.markers_summary.tsv.json') as j:
            assert json.load(j) == {'data': [
                {'user_genome': 'somefile1.fasta', 'field1': 'fee', 'field2': 'fie'},
                {'user_genome': 'somefile2.fasta', 'field1': 'fo', 'field2': 'fum'},
            ]}
