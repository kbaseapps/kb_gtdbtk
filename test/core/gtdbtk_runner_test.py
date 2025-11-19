import json
import logging
import os
import tempfile
from unittest.mock import Mock
import pytest

from pathlib import Path

from kb_gtdbtk.core.gtdbtk_runner import (
    run_gtdbtk,
    _load_summary_tsv_file,
    _merge_summary_tsv_files,
    _process_output_files
)

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s', level=logging.INFO)


def test_gtdbtk_run():
    db_ver = 214

    with tempfile.TemporaryDirectory(prefix='test_gtdbtk_run') as test_dir_str:
        test_dir = Path(test_dir_str)
        out_dir = test_dir / 'output'
        out_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = test_dir / 'temp'
        temp_dir.mkdir(parents=True, exist_ok=True)
        data_dir = test_dir / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
        db_dir = data_dir / f"r{db_ver}" / "mash"
        db_dir.mkdir(parents=True, exist_ok=True)
        refdata_file = db_dir / 'gtdb_ref_sketch.msh'
        refdata_file.touch(exist_ok=True)

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

            expected_cmd = [
                "gtdbtk",
                "classify_wf",
                "--out_dir", str(temp_out),
                "--batchfile",
                "--cpus", "16",
                "--min_perc_aa", "50.2"
            ]

            if "--skip_ani_screen" in command:
                expected_cmd += [
                    "--skip_ani_screen",
                    "--no_mash"
                ]
            else:
                expected_cmd += [
                    "--mash_db",
                    str(data_dir / f"r{db_ver}" / "mash" / "gtdb_ref_sketch.msh")
                ]

            assert command == expected_cmd

            # arbitrary TSV files, these do not match what GTDB-tk produces
            # summary must have the correct number of fields

            with open(os.path.join(temp_classify, 'gtdbtk.ar53.summary.tsv'), 'w') as t:
                t.writelines(['\t'.join(['user_genome', 'classification', 'fastani_reference', 'fastani_reference_radius', 'fastani_taxonomy', 'fastani_ani', 'fastani_af', 'closest_placement_reference', 'closest_placement_radius', 'closest_placement_taxonomy', 'closest_placement_ani', 'closest_placement_af', 'pplacer_taxonomy', 'classification_method', 'note', 'other_related_references', 'msa_percent', 'translation_table', 'red_value', 'warnings']) + '\n',
                              '\t'.join(['id0', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo']) + '\n',
                              '\t'.join(['id1', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo']) + '\n'
                              ])

            with open(os.path.join(temp_classify, 'gtdbtk.bac120.summary.tsv'), 'w') as t:
                t.writelines(['\t'.join(['user_genome', 'classification', 'fastani_reference', 'fastani_reference_radius', 'fastani_taxonomy', 'fastani_ani', 'fastani_af', 'closest_placement_reference', 'closest_placement_radius', 'closest_placement_taxonomy', 'closest_placement_ani', 'closest_placement_af', 'pplacer_taxonomy', 'classification_method', 'note', 'other_related_references', 'msa_percent', 'translation_table', 'red_value', 'warnings']) + '\n',
                              '\t'.join(['id0', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo']) + '\n',
                              '\t'.join(['id1', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo']) + '\n'
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
            {
                Path('/somepath1'): 'somefile1.fasta',
                Path('/somepath2'): 'somefile2.fasta',
            },
            out_dir,
            temp_dir,
            50.2,
            db_ver,
            0,
            16,
            data_root_dir=data_dir
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
                {'user_genome': 'somefile1.fasta', 'classification':'foo', 'fastani_reference':'foo', 'fastani_reference_radius':'foo', 'fastani_taxonomy':'foo', 'fastani_ani':'foo', 'fastani_af':'foo', 'closest_placement_reference':'foo', 'closest_placement_radius':'foo', 'closest_placement_taxonomy':'foo', 'closest_placement_ani':'foo', 'closest_placement_af':'foo', 'pplacer_taxonomy':'foo', 'classification_method':'foo', 'note':'foo', 'other_related_references':'foo', 'msa_percent':'foo', 'translation_table':'foo', 'red_value':'foo', 'warnings':'foo'},
                {'user_genome': 'somefile2.fasta', 'classification':'foo', 'fastani_reference':'foo', 'fastani_reference_radius':'foo', 'fastani_taxonomy':'foo', 'fastani_ani':'foo', 'fastani_af':'foo', 'closest_placement_reference':'foo', 'closest_placement_radius':'foo', 'closest_placement_taxonomy':'foo', 'closest_placement_ani':'foo', 'closest_placement_af':'foo', 'pplacer_taxonomy':'foo', 'classification_method':'foo', 'note':'foo', 'other_related_references':'foo', 'msa_percent':'foo', 'translation_table':'foo', 'red_value':'foo', 'warnings':'foo'}
            ]}

        with open(os.path.join(out_dir, 'gtdbtk.bac120.summary.tsv.json')) as j:
            assert json.load(j) == {'data': [
                {'user_genome': 'somefile1.fasta', 'classification':'foo', 'fastani_reference':'foo', 'fastani_reference_radius':'foo', 'fastani_taxonomy':'foo', 'fastani_ani':'foo', 'fastani_af':'foo', 'closest_placement_reference':'foo', 'closest_placement_radius':'foo', 'closest_placement_taxonomy':'foo', 'closest_placement_ani':'foo', 'closest_placement_af':'foo', 'pplacer_taxonomy':'foo', 'classification_method':'foo', 'note':'foo', 'other_related_references':'foo', 'msa_percent':'foo', 'translation_table':'foo', 'red_value':'foo', 'warnings':'foo'},
                {'user_genome': 'somefile2.fasta', 'classification':'foo', 'fastani_reference':'foo', 'fastani_reference_radius':'foo', 'fastani_taxonomy':'foo', 'fastani_ani':'foo', 'fastani_af':'foo', 'closest_placement_reference':'foo', 'closest_placement_radius':'foo', 'closest_placement_taxonomy':'foo', 'closest_placement_ani':'foo', 'closest_placement_af':'foo', 'pplacer_taxonomy':'foo', 'classification_method':'foo', 'note':'foo', 'other_related_references':'foo', 'msa_percent':'foo', 'translation_table':'foo', 'red_value':'foo', 'warnings':'foo'}
            ]}

        with open(os.path.join(out_dir, 'gtdbtk.bac120.markers_summary.tsv.json')) as j:
            assert json.load(j) == {'data': [
                {'user_genome': 'somefile1.fasta', 'field1': 'fee', 'field2': 'fie'},
                {'user_genome': 'somefile2.fasta', 'field1': 'fo', 'field2': 'fum'},
            ]}

def test_gtdbtk_run_fail_no_default_refdata(tmp_path):
    db_ver = 214
    expected_db_path = f"/data/r{db_ver}/mash/gtdb_ref_sketch.msh"
    with pytest.raises(RuntimeError, match=f"GTDB ref genomes MASH DB not found in expected path {expected_db_path}. This must be generated during refdata initialization."):
        run_gtdbtk(
            Mock(),
            {
                Path('/somepath1'): 'somefile1.fasta',
                Path('/somepath2'): 'somefile2.fasta',
            },
            tmp_path / "out_dir",
            tmp_path / "temp_dir",
            50.2,
            db_ver,
            0,
            16,
        )


class TestLoadSummaryTsvFile:
    """Tests for _load_summary_tsv_file function"""

    def test_load_summary_tsv_basic(self, tmp_path):
        """Test loading a basic TSV file with header and data rows"""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("user_genome\tfield1\tfield2\n"
                           "genome1\tvalue1\tvalue2\n"
                           "genome2\tvalue3\tvalue4\n")

        result = _load_summary_tsv_file(tsv_file)

        assert result["header"] == ["user_genome", "field1", "field2"]
        assert result["data"]["genome1"] == ["genome1", "value1", "value2"]
        assert result["data"]["genome2"] == ["genome2", "value3", "value4"]
        assert result["id_order"] == ["genome1", "genome2"]

    def test_load_summary_tsv_with_na_values(self, tmp_path):
        """Test loading TSV with N/A values"""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("user_genome\tfield1\tfield2\n"
                           "genome1\tN/A\tvalue2\n"
                           "genome2\tvalue3\tN/A\n")

        result = _load_summary_tsv_file(tsv_file)

        assert result["data"]["genome1"][1] == "N/A"
        assert result["data"]["genome2"][2] == "N/A"

    def test_load_summary_tsv_preserves_trailing_tabs(self, tmp_path):
        """Test that trailing empty fields are preserved"""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("user_genome\tfield1\tfield2\tfield3\n"
                           "genome1\tvalue1\tvalue2\t\n")

        result = _load_summary_tsv_file(tsv_file)

        # The last element should be empty string, not removed
        assert len(result["data"]["genome1"]) == 4
        assert result["data"]["genome1"][3] == ""

class TestMergeSummaryTsvFiles:
    """Tests for _merge_summary_tsv_files function"""

    def test_merge_both_files_exist(self, tmp_path):
        """Test merging when both std and tree files exist"""
        std_file = tmp_path / "std.tsv"
        tree_file = tmp_path / "tree.tsv"
        out_file = tmp_path / "merged.tsv"

        # Standard file with some N/A values
        std_file.write_text("user_genome\tfield1\tfield2\tfield3\n"
                           "genome1\tN/A\tvalue2\tvalue3\n"
                           "genome2\tvalue1\tN/A\tvalue3\n")

        # Tree file with values where std has N/A
        tree_file.write_text("user_genome\tfield1\tfield2\tfield3\n"
                            "genome1\tvalue1_tree\tN/A\tvalue3\n"
                            "genome2\tvalue1_tree\tvalue2_tree\tN/A\n")

        _merge_summary_tsv_files(std_file, tree_file, out_file)

        result = _load_summary_tsv_file(out_file)

        # Tree values should override N/A in std
        assert result["data"]["genome1"][1] == "value1_tree"  # Was N/A in std
        assert result["data"]["genome1"][2] == "value2"  # From std, tree has N/A
        assert result["data"]["genome2"][1] == "value1"  # From std
        assert result["data"]["genome2"][2] == "value2_tree"  # Was N/A in std

    def test_merge_only_std_file(self, tmp_path):
        """Test when only standard file exists"""
        std_file = tmp_path / "std.tsv"
        tree_file = tmp_path / "tree.tsv"
        out_file = tmp_path / "merged.tsv"

        std_file.write_text("user_genome\tfield1\tfield2\n"
                           "genome1\tvalue1\tvalue2\n")

        _merge_summary_tsv_files(std_file, tree_file, out_file)

        assert out_file.read_text() == std_file.read_text()

    def test_merge_only_tree_file(self, tmp_path):
        """Test when only tree file exists"""
        std_file = tmp_path / "std.tsv"
        tree_file = tmp_path / "tree.tsv"
        out_file = tmp_path / "merged.tsv"

        tree_file.write_text("user_genome\tfield1\tfield2\n"
                            "genome1\tvalue1\tvalue2\n")

        _merge_summary_tsv_files(std_file, tree_file, out_file)

        assert out_file.read_text() == tree_file.read_text()

    def test_merge_no_files(self, tmp_path):
        """Test when neither file exists"""
        std_file = tmp_path / "std.tsv"
        tree_file = tmp_path / "tree.tsv"
        out_file = tmp_path / "merged.tsv"

        _merge_summary_tsv_files(std_file, tree_file, out_file)

        # Output file should not be created
        assert not out_file.exists()

    def test_merge_mismatched_headers(self, tmp_path):
        """Test merging files with different headers"""
        std_file = tmp_path / "std.tsv"
        tree_file = tmp_path / "tree.tsv"
        out_file = tmp_path / "merged.tsv"

        std_file.write_text("user_genome\tfield1\tfield2\n"
                           "genome1\tvalue1\tvalue2\n")

        tree_file.write_text("user_genome\tfield1\tfield3\n"
                            "genome1\tvalue1_tree\tvalue3\n")

        _merge_summary_tsv_files(std_file, tree_file, out_file)

        # Should copy only the tree file when headers don't match
        assert out_file.read_text() == tree_file.read_text()


class TestProcessOutputFiles:
    """Tests for _process_output_files function"""

    def test_process_output_creates_json_files(self, tmp_path):
        """Test that JSON versions of TSV files are created"""
        temp_output = tmp_path / "temp_output"
        temp_output.mkdir()
        (temp_output / "classify").mkdir()
        (temp_output / "identify").mkdir()

        temp_trees_output = tmp_path / "temp_trees_output"
        temp_trees_output.mkdir()

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        # Create a summary TSV file
        summary_file = temp_output / "classify" / "gtdbtk.bac120.summary.tsv"
        summary_file.write_text("user_genome\tclassification\tfield2\n"
                               "id0\tclass_a\tvalue1\n"
                               "id1\tclass_b\tvalue2\n")

        id_to_name = {
            "id0": "genome1.fasta",
            "id1": "genome2.fasta"
        }

        classification, summary_tables = _process_output_files(
            temp_output, temp_trees_output, out_dir, id_to_name
        )

        # Check that JSON file was created
        json_file = out_dir / "gtdbtk.bac120.summary.tsv.json"
        assert json_file.exists()

        # Check JSON content
        with open(json_file) as f:
            json_data = json.load(f)

        assert "data" in json_data
        assert len(json_data["data"]) == 2
        assert json_data["data"][0]["user_genome"] == "genome1.fasta"
        assert json_data["data"][1]["user_genome"] == "genome2.fasta"

    def test_process_output_classification_mapping(self, tmp_path):
        """Test that classification dict is populated correctly"""
        temp_output = tmp_path / "temp_output"
        temp_output.mkdir()
        (temp_output / "classify").mkdir()
        (temp_output / "identify").mkdir()

        temp_trees_output = tmp_path / "temp_trees_output"
        temp_trees_output.mkdir()

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        summary_file = temp_output / "classify" / "gtdbtk.bac120.summary.tsv"
        summary_file.write_text("user_genome\tclassification\n"
                               "id0\td__Bacteria;p__Test\n"
                               "id1\td__Bacteria;p__Other\n")

        id_to_name = {
            "id0": "genome1.fasta",
            "id1": "genome2.fasta"
        }

        classification, _ = _process_output_files(
            temp_output, temp_trees_output, out_dir, id_to_name
        )

        assert classification["genome1.fasta"] == "d__Bacteria;p__Test"
        assert classification["genome2.fasta"] == "d__Bacteria;p__Other"

    def test_process_output_handles_null_values(self, tmp_path):
        """Test that null/blank values are replaced with '-'"""
        temp_output = tmp_path / "temp_output"
        temp_output.mkdir()
        (temp_output / "classify").mkdir()
        (temp_output / "identify").mkdir()

        temp_trees_output = tmp_path / "temp_trees_output"
        temp_trees_output.mkdir()

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        summary_file = temp_output / "classify" / "gtdbtk.bac120.summary.tsv"
        summary_file.write_text("user_genome\tfield1\tfield2\n"
                               "id0\tvalue1\t\n"
                               "id1\t\tvalue2\n")

        id_to_name = {
            "id0": "genome1.fasta",
            "id1": "genome2.fasta"
        }

        _, summary_tables = _process_output_files(
            temp_output, temp_trees_output, out_dir, id_to_name
        )

        json_data = summary_tables["gtdbtk.bac120.summary.tsv"]
        # Empty values should be replaced with '-'
        assert json_data["data"][0]["field2"] == "-"
        assert json_data["data"][1]["field1"] == "-"

    def test_process_output_creates_id_map_file(self, tmp_path):
        """Test that id_to_name.map file is created"""
        temp_output = tmp_path / "temp_output"
        temp_output.mkdir()
        (temp_output / "classify").mkdir()
        (temp_output / "identify").mkdir()

        temp_trees_output = tmp_path / "temp_trees_output"
        temp_trees_output.mkdir()

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        id_to_name = {
            "id0": "genome1.fasta",
            "id1": "genome2.fasta"
        }

        _process_output_files(
            temp_output, temp_trees_output, out_dir, id_to_name
        )

        map_file = out_dir / "id_to_name.map"
        assert map_file.exists()

        content = map_file.read_text()
        assert "id0\tgenome1.fasta\n" in content
        assert "id1\tgenome2.fasta\n" in content

    def test_process_output_missing_id_in_mapping_raises_error(self, tmp_path):
        """Test that missing ID in mapping raises ValueError"""
        temp_output = tmp_path / "temp_output"
        temp_output.mkdir()
        (temp_output / "classify").mkdir()
        (temp_output / "identify").mkdir()

        temp_trees_output = tmp_path / "temp_trees_output"
        temp_trees_output.mkdir()

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        summary_file = temp_output / "classify" / "gtdbtk.bac120.summary.tsv"
        summary_file.write_text("user_genome\tclassification\n"
                               "id_missing\tclass_a\n")

        id_to_name = {
            "id0": "genome1.fasta"
        }

        with pytest.raises(ValueError, match="missing id_missing"):
            _process_output_files(
                temp_output, temp_trees_output, out_dir, id_to_name
            )
