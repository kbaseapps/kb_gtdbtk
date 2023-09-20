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

            #temp_out = temp_dir / 'output'
            temp_out = str(td)
            temp_classify = Path(os.path.join(temp_out, 'classify'))
            temp_classify.mkdir(parents=True, exist_ok=True)
            temp_identify = Path(os.path.join(temp_out, 'identify'))
            temp_identify.mkdir(parents=True, exist_ok=True)
            
            assert command == [
                'gtdbtk',
                'classify_wf',
                '--out_dir', str(temp_out),
                '--batchfile',  # arg popped
                '--cpus', '16',
                '--min_perc_aa', '50.2',
                '--skip_ani_screen']

            # arbitrary TSV files, these may not match what GTDB-tk produces
            # implmentation doesn't care for now other than the first column

            with open(os.path.join(temp_classify, 'gtdbtk.ar53.summary.tsv'), 'w') as t:
                t.writelines(['\t'.join(['name', 'field1', 'field2']) + '\n',
                              '\t'.join(['id0', 'foo', 'bar']) + '\n',
                              '\t'.join(['id1', 'baz', 'bat']) + '\n',
                              ])

            with open(os.path.join(temp_classify, 'gtdbtk.bac120.summary.tsv'), 'w') as t:
                t.writelines(['\t'.join(['name', 'field1', 'field2']) + '\n',
                              '\t'.join(['id0', 'whoo', 'whee']) + '\n',
                              '\t'.join(['id1', 'whoa', 'whump']) + '\n',
                              ])

            with open(os.path.join(temp_identify, 'gtdbtk.bac120.markers_summary.tsv'), 'w') as t:
                t.writelines(['\t'.join(['user_genome', 'field1', 'field2']) + '\n',
                              '\t'.join(['id0', 'fee', 'fie']) + '\n',
                              '\t'.join(['id1', 'fo', 'fum']) + '\n',
                              ])

            # skip 'gtdbtk.ar53.markers_summary.tsv'

            # ####   end of runner callable   ####

        run_gtdbtk(
            runner,
            {Path('/somepath1'): 'somefile1.fasta',
             Path('/somepath2'): 'somefile2.fasta',
             },
            out_dir,
            temp_dir,
            50.2,
            207,
            0,
            16
            )

        with open(tf[0]) as bf:
            lines = bf.readlines()
            assert len(lines) == 2
            #assert lines[0] == f'{temp_dir}/links/id0\tid0\n'  # made path unique with timestamp
            #assert lines[1] == f'{temp_dir}/links/id1\tid1\n'
            id0_path = lines[0].split("\t")[0]
            id1_path = lines[1].split("\t")[0]
            
        assert os.readlink(id0_path) == '/somepath1'
        assert os.readlink(id1_path) == '/somepath2'

        assert sorted(os.listdir(out_dir)) == [
            #'gtdbtk.ar53.classify.tree',
            #'gtdbtk.ar53.markers_summary.tsv',
            #'gtdbtk.ar53.markers_summary.tsv.json',
            'gtdbtk.ar53.summary.tsv',
            'gtdbtk.ar53.summary.tsv.json',
            #'gtdbtk.bac120.classify.tree',
            'gtdbtk.bac120.markers_summary.tsv',
            'gtdbtk.bac120.markers_summary.tsv.json',
            'gtdbtk.bac120.summary.tsv',
            'gtdbtk.bac120.summary.tsv.json',
            'id_to_name.map',
            'runtime_output'
        ]

        with open(os.path.join(out_dir, 'gtdbtk.ar53.summary.tsv.json')) as j:
            assert json.load(j) == {'data': [
                {'name': 'somefile1.fasta', 'field1': 'foo', 'field2': 'bar'},
                {'name': 'somefile2.fasta', 'field1': 'baz', 'field2': 'bat'},
            ]}

        with open(os.path.join(out_dir, 'gtdbtk.bac120.summary.tsv.json')) as j:
            assert json.load(j) == {'data': [
                {'name': 'somefile1.fasta', 'field1': 'whoo', 'field2': 'whee'},
                {'name': 'somefile2.fasta', 'field1': 'whoa', 'field2': 'whump'},
            ]}

        with open(os.path.join(out_dir, 'gtdbtk.bac120.markers_summary.tsv.json')) as j:
            assert json.load(j) == {'data': [
                {'user_genome': 'somefile1.fasta', 'field1': 'fee', 'field2': 'fie'},
                {'user_genome': 'somefile2.fasta', 'field1': 'fo', 'field2': 'fum'},
            ]}
