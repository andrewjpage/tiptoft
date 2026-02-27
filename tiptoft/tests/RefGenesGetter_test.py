"""
Tests for tiptoft/RefGenesGetter.py

Covers:
  - Error exception class
  - RefGenesGetter initialisation (default and verbose)
  - Regression guard: corrected filename 'enterobacteriales.fsa' is used
  - _validate_fasta_file: valid FASTA passes, HTML/empty/missing/edge-case files
  - run() dispatches correctly and forwards outprefix
  - _get_from_plasmidfinder integration (subprocess, validate, shutil mocked)
  - Duplicate sequence ID de-duplication logic
  - ID normalisation (underscore → dot replacement)
  - Sequence uppercasing
  - TSV metadata output format
  - URL construction and base URL correctness
  - Non-.fsa files ignored during combining
  - subprocess failure propagation
"""

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from tiptoft.RefGenesGetter import Error, RefGenesGetter


# ---------------------------------------------------------------------------
# TestErrorClass
# ---------------------------------------------------------------------------

class TestErrorClass(unittest.TestCase):
    """Tests for the custom Error exception class."""

    def test_error_is_exception_subclass(self):
        self.assertTrue(issubclass(Error, Exception))

    def test_error_message_is_preserved(self):
        err = Error('something went wrong')
        self.assertEqual(str(err), 'something went wrong')

    def test_error_can_be_raised_and_caught(self):
        with self.assertRaises(Error) as ctx:
            raise Error('test message')
        self.assertEqual(str(ctx.exception), 'test message')


# ---------------------------------------------------------------------------
# TestRefGenesGetterInit
# ---------------------------------------------------------------------------

class TestRefGenesGetterInit(unittest.TestCase):
    """Tests for __init__ default and verbose modes."""

    def test_default_initialization(self):
        getter = RefGenesGetter()
        self.assertEqual(getter.ref_db, 'plasmidfinder')
        self.assertFalse(getter.verbose)
        self.assertEqual(getter.genetic_code, 11)
        self.assertEqual(getter.max_download_attempts, 3)
        self.assertEqual(getter.sleep_time, 2)

    def test_verbose_initialization(self):
        getter = RefGenesGetter(verbose=True)
        self.assertTrue(getter.verbose)


# ---------------------------------------------------------------------------
# TestFilesToDownload
# ---------------------------------------------------------------------------

class TestFilesToDownload(unittest.TestCase):
    """Regression guard: correct filenames are passed to curl."""

    EXPECTED_FILES = [
        'enterobacteriales.fsa',
        'Inc18.fsa',
        'NT_Rep.fsa',
        'Rep1.fsa',
        'Rep2.fsa',
        'Rep3.fsa',
        'RepA_N.fsa',
        'RepL.fsa',
        'Rep_trans.fsa',
    ]

    def _collect_downloaded_filenames(self):
        """Run _get_from_plasmidfinder with mocked I/O and return filenames passed to curl."""
        downloaded = []

        def fake_check_call(cmd, shell=False):
            # cmd format: 'curl -o <filename> <url>'
            parts = cmd.split()
            downloaded.append(parts[2])

        getter = RefGenesGetter()
        with patch('subprocess.check_call', side_effect=fake_check_call), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[]), \
             patch('os.listdir', return_value=[]), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        return downloaded

    def test_enterobacteriales_is_in_download_list(self):
        downloaded = self._collect_downloaded_filenames()
        self.assertIn('enterobacteriales.fsa', downloaded)

    def test_enterobacteriaceae_is_not_in_download_list(self):
        downloaded = self._collect_downloaded_filenames()
        self.assertNotIn('enterobacteriaceae.fsa', downloaded)

    def test_all_expected_filenames_are_present(self):
        downloaded = self._collect_downloaded_filenames()
        self.assertEqual(sorted(downloaded), sorted(self.EXPECTED_FILES))

    def test_download_count_is_correct(self):
        downloaded = self._collect_downloaded_filenames()
        self.assertEqual(len(downloaded), 9)


# ---------------------------------------------------------------------------
# TestValidateFastaFile
# ---------------------------------------------------------------------------

class TestValidateFastaFile(unittest.TestCase):
    """Tests for _validate_fasta_file using real temp files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.getter = RefGenesGetter()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_file(self, content, filename='test.fsa'):
        path = os.path.join(self.tmpdir, filename)
        with open(path, 'w') as fh:
            fh.write(content)
        return path

    def test_valid_fasta_passes_validation(self):
        path = self._make_file('>seq1\nATGC\n')
        self.getter._validate_fasta_file(path)  # should not raise

    def test_html_doctype_raises_error(self):
        path = self._make_file('<!DOCTYPE html>\n<html></html>\n')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(path)

    def test_html_tag_raises_error(self):
        path = self._make_file('<html>\n<body>Not found</body>\n</html>\n')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(path)

    def test_empty_file_raises_error(self):
        path = self._make_file('')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(path)

    def test_nonexistent_file_raises_error(self):
        missing = os.path.join(self.tmpdir, 'does_not_exist.fsa')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(missing)

    def test_error_message_includes_filepath(self):
        path = self._make_file('<!DOCTYPE html>\n')
        with self.assertRaises(Error) as ctx:
            self.getter._validate_fasta_file(path)
        self.assertIn(path, str(ctx.exception))

    def test_error_message_includes_first_line(self):
        path = self._make_file('<!DOCTYPE html>\n<html>\n')
        with self.assertRaises(Error) as ctx:
            self.getter._validate_fasta_file(path)
        self.assertIn('<!DOCTYPE html>', str(ctx.exception))

    def test_xml_declaration_raises_error(self):
        """XML responses also start with '<' and should be caught."""
        path = self._make_file('<?xml version="1.0"?>\n<error/>\n')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(path)

    def test_html_with_leading_whitespace_raises_error(self):
        """Leading whitespace before '<' should still be detected after strip()."""
        path = self._make_file('  <!DOCTYPE html>\n')
        with self.assertRaises(Error):
            self.getter._validate_fasta_file(path)

    def test_newline_only_file_passes_validation(self):
        """A file with only a newline is not empty (readline returns '\\n'), not HTML."""
        path = self._make_file('\n')
        # first_line = '\n', which is truthy, and '\n'.strip() = '' which doesn't start with '<'
        self.getter._validate_fasta_file(path)  # should not raise

    def test_fasta_with_description_passes(self):
        """Multi-word FASTA header lines are valid."""
        path = self._make_file('>seq1 some description here\nATGCATGC\n')
        self.getter._validate_fasta_file(path)  # should not raise

    def test_empty_error_mentions_empty(self):
        """Empty file error message should contain 'empty'."""
        path = self._make_file('')
        with self.assertRaises(Error) as ctx:
            self.getter._validate_fasta_file(path)
        self.assertIn('empty', str(ctx.exception).lower())

    def test_nonexistent_error_mentions_could_not_read(self):
        """Missing file error message should mention read failure."""
        missing = os.path.join(self.tmpdir, 'nope.fsa')
        with self.assertRaises(Error) as ctx:
            self.getter._validate_fasta_file(missing)
        self.assertIn('Could not read', str(ctx.exception))


# ---------------------------------------------------------------------------
# TestRunMethod
# ---------------------------------------------------------------------------

class TestRunMethod(unittest.TestCase):
    """Tests that run() dispatches to _get_from_plasmidfinder correctly."""

    def test_run_dispatches_to_plasmidfinder(self):
        getter = RefGenesGetter()
        with patch.object(getter, '_get_from_plasmidfinder') as mock_method:
            getter.run('some_prefix')
        mock_method.assert_called_once()

    def test_run_passes_outprefix_unchanged(self):
        getter = RefGenesGetter()
        with patch.object(getter, '_get_from_plasmidfinder') as mock_method:
            getter.run('my_output_prefix')
        mock_method.assert_called_once_with('my_output_prefix')


# ---------------------------------------------------------------------------
# TestGetFromPlasmidfinderIntegration
# ---------------------------------------------------------------------------

class TestGetFromPlasmidfinderIntegration(unittest.TestCase):
    """Integration tests for _get_from_plasmidfinder with I/O mocked."""

    def _run_with_mocks(self, verbose=False, validate_side_effect=None):
        """Helper: run _get_from_plasmidfinder with all external I/O mocked."""
        getter = RefGenesGetter(verbose=verbose)

        mock_validate = MagicMock(side_effect=validate_side_effect)
        mock_check_call = MagicMock()
        mock_open_write = MagicMock(return_value=MagicMock())
        mock_close = MagicMock()
        mock_file_reader = MagicMock(return_value=[])
        mock_listdir = MagicMock(return_value=[])
        mock_rmtree = MagicMock()

        with patch('subprocess.check_call', mock_check_call), \
             patch.object(getter, '_validate_fasta_file', mock_validate), \
             patch('pyfastaq.utils.open_file_write', mock_open_write), \
             patch('pyfastaq.utils.close', mock_close), \
             patch('pyfastaq.sequences.file_reader', mock_file_reader), \
             patch('os.listdir', mock_listdir), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('shutil.rmtree', mock_rmtree):
            getter._get_from_plasmidfinder('test_out')

        return getter, mock_check_call, mock_validate, mock_rmtree

    def test_curl_is_called_for_each_file(self):
        _, mock_check_call, _, _ = self._run_with_mocks()
        self.assertEqual(mock_check_call.call_count, 9)

    def test_validate_is_called_for_each_downloaded_file(self):
        _, _, mock_validate, _ = self._run_with_mocks()
        self.assertEqual(mock_validate.call_count, 9)

    def test_tmpdir_removed_when_not_verbose(self):
        getter, _, _, mock_rmtree = self._run_with_mocks(verbose=False)
        outprefix = os.path.abspath('test_out')
        expected_tmpdir = outprefix + '.tmp.download'
        mock_rmtree.assert_called_once_with(expected_tmpdir)

    def test_tmpdir_retained_when_verbose(self):
        _, _, _, mock_rmtree = self._run_with_mocks(verbose=True)
        mock_rmtree.assert_not_called()

    def test_html_response_propagates_error(self):
        side_effect = Error('Downloaded file appears to be an HTML error page')
        with self.assertRaises(Error):
            self._run_with_mocks(validate_side_effect=side_effect)

    def test_curl_url_contains_correct_filename(self):
        _, mock_check_call, _, _ = self._run_with_mocks()
        curl_cmds = [args[0][0] for args in mock_check_call.call_args_list]

        enterobacteriales_cmds = [c for c in curl_cmds if 'enterobacteriales' in c]
        enterobacteriaceae_cmds = [c for c in curl_cmds if 'enterobacteriaceae' in c]
        self.assertEqual(len(enterobacteriales_cmds), 1)
        self.assertEqual(len(enterobacteriaceae_cmds), 0)

    def test_curl_urls_use_correct_base_url(self):
        """All curl commands should target the PlasmidFinder BitBucket repo."""
        _, mock_check_call, _, _ = self._run_with_mocks()
        expected_base = 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/master/'
        for call_args in mock_check_call.call_args_list:
            cmd = call_args[0][0]
            self.assertIn(expected_base, cmd)

    def test_subprocess_failure_propagates(self):
        """CalledProcessError from curl should propagate uncaught."""
        getter = RefGenesGetter()
        mock_check_call = MagicMock(
            side_effect=subprocess.CalledProcessError(1, 'curl'))

        with patch('subprocess.check_call', mock_check_call), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[]), \
             patch('os.listdir', return_value=[]), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('shutil.rmtree'):
            with self.assertRaises(subprocess.CalledProcessError):
                getter._get_from_plasmidfinder('test_out')

    def test_non_fsa_files_are_ignored(self):
        """Files without .fsa extension in tmpdir should not be processed."""
        getter = RefGenesGetter()
        mock_file_reader = MagicMock(return_value=[])

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', mock_file_reader), \
             patch('os.listdir', return_value=['readme.txt', 'data.csv', 'real.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        # file_reader should only be called for the .fsa file
        self.assertEqual(mock_file_reader.call_count, 1)

    def test_outprefix_is_converted_to_absolute_path(self):
        """Output files should use absolute paths even when relative prefix given."""
        getter = RefGenesGetter()
        mock_open_write = MagicMock(return_value=MagicMock())

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', mock_open_write), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[]), \
             patch('os.listdir', return_value=[]), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('relative_prefix')

        # Both .fa and .tsv open_file_write calls should use absolute paths
        for call_args in mock_open_write.call_args_list:
            path = call_args[0][0]
            self.assertTrue(os.path.isabs(path),
                            f'Expected absolute path, got: {path}')


# ---------------------------------------------------------------------------
# TestDuplicateIdHandling
# ---------------------------------------------------------------------------

class TestDuplicateIdHandling(unittest.TestCase):
    """Tests that duplicate sequence IDs are correctly de-duplicated."""

    def _run_with_sequences(self, sequences):
        """Helper: run _get_from_plasmidfinder with given sequences and capture written IDs."""
        getter = RefGenesGetter()
        written_ids = []

        def fake_print(*args, **kwargs):
            # Capture sequence objects written via `print(seq, file=fout_fa)`
            if len(args) == 1 and not isinstance(args[0], str):
                try:
                    written_ids.append(args[0].id)
                except AttributeError:
                    pass

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=sequences), \
             patch('os.listdir', return_value=['test.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('builtins.print', side_effect=fake_print), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        return written_ids

    def _make_seq(self, seq_id, seq_data='ATGC'):
        seq = MagicMock()
        seq.id = seq_id
        seq.seq = seq_data
        return seq

    def test_duplicate_ids_are_disambiguated(self):
        """Two sequences with 'rep_1' become 'rep.1' and 'rep.1.2'."""
        ids = self._run_with_sequences([
            self._make_seq('rep_1'),
            self._make_seq('rep_1'),
        ])
        self.assertIn('rep.1', ids)
        self.assertIn('rep.1.2', ids)

    def test_three_duplicates_get_incremented_suffixes(self):
        """Three identical IDs get .1 (implicit), .1.2, .1.3."""
        ids = self._run_with_sequences([
            self._make_seq('abc_x'),
            self._make_seq('abc_x'),
            self._make_seq('abc_x'),
        ])
        self.assertIn('abc.x', ids)
        self.assertIn('abc.x.2', ids)
        self.assertIn('abc.x.3', ids)

    def test_unique_ids_are_not_modified_beyond_normalization(self):
        """Unique IDs only get underscore→dot normalization, no suffix."""
        ids = self._run_with_sequences([
            self._make_seq('alpha_1'),
            self._make_seq('beta_2'),
        ])
        self.assertIn('alpha.1', ids)
        self.assertIn('beta.2', ids)
        self.assertEqual(len(ids), 2)


# ---------------------------------------------------------------------------
# TestIdNormalization
# ---------------------------------------------------------------------------

class TestIdNormalization(unittest.TestCase):
    """Tests for the ID normalization logic (underscore → dot replacement)."""

    def _run_with_sequences(self, sequences):
        """Helper: run _get_from_plasmidfinder and capture written IDs."""
        getter = RefGenesGetter()
        written_ids = []

        def fake_print(*args, **kwargs):
            if len(args) == 1 and not isinstance(args[0], str):
                try:
                    written_ids.append(args[0].id)
                except AttributeError:
                    pass

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=sequences), \
             patch('os.listdir', return_value=['test.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('builtins.print', side_effect=fake_print), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        return written_ids

    def _make_seq(self, seq_id, seq_data='ATGC'):
        seq = MagicMock()
        seq.id = seq_id
        seq.seq = seq_data
        return seq

    def test_first_underscore_replaced_with_dot(self):
        """'rep_1' → 'rep.1'."""
        ids = self._run_with_sequences([self._make_seq('rep_1')])
        self.assertEqual(ids, ['rep.1'])

    def test_only_first_underscore_is_replaced(self):
        """'rep_1_extra' → 'rep.1_extra' (only first underscore)."""
        ids = self._run_with_sequences([self._make_seq('rep_1_extra')])
        self.assertEqual(ids, ['rep.1_extra'])

    def test_no_underscore_id_unchanged(self):
        """An ID without underscores passes through unchanged."""
        ids = self._run_with_sequences([self._make_seq('plasmid123')])
        self.assertEqual(ids, ['plasmid123'])

    def test_sequence_is_uppercased(self):
        """Sequences should be uppercased during processing."""
        seq = self._make_seq('test_1', 'atgc')
        getter = RefGenesGetter()

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[seq]), \
             patch('os.listdir', return_value=['test.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('builtins.print', lambda *a, **kw: None), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        self.assertEqual(seq.seq, 'ATGC')


# ---------------------------------------------------------------------------
# TestTsvOutput
# ---------------------------------------------------------------------------

class TestTsvOutput(unittest.TestCase):
    """Tests for TSV metadata output format."""

    def test_tsv_line_contains_original_id(self):
        """TSV output should include 'Original name was <id>' with the pre-normalized ID."""
        getter = RefGenesGetter()
        seq = MagicMock()
        seq.id = 'rep_1'
        seq.seq = 'ATGC'

        tsv_lines = []

        def fake_print(*args, **kwargs):
            file_arg = kwargs.get('file')
            sep = kwargs.get('sep', ' ')
            # TSV prints use sep='\t' and have 6+ positional args
            if sep == '\t' and len(args) >= 6:
                tsv_lines.append(sep.join(str(a) for a in args))

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[seq]), \
             patch('os.listdir', return_value=['test.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('builtins.print', side_effect=fake_print), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        self.assertEqual(len(tsv_lines), 1)
        self.assertIn('Original name was rep_1', tsv_lines[0])

    def test_tsv_has_six_tab_separated_columns(self):
        """Each TSV line should have exactly 6 tab-separated fields."""
        getter = RefGenesGetter()
        seq = MagicMock()
        seq.id = 'test_seq'
        seq.seq = 'ATGC'

        tsv_lines = []

        def fake_print(*args, **kwargs):
            sep = kwargs.get('sep', ' ')
            if sep == '\t' and len(args) >= 6:
                tsv_lines.append(sep.join(str(a) for a in args))

        with patch('subprocess.check_call'), \
             patch.object(getter, '_validate_fasta_file'), \
             patch('pyfastaq.utils.open_file_write', return_value=MagicMock()), \
             patch('pyfastaq.utils.close'), \
             patch('pyfastaq.sequences.file_reader', return_value=[seq]), \
             patch('os.listdir', return_value=['test.fsa']), \
             patch('os.mkdir'), \
             patch('os.chdir'), \
             patch('builtins.print', side_effect=fake_print), \
             patch('shutil.rmtree'):
            getter._get_from_plasmidfinder('test_out')

        self.assertEqual(len(tsv_lines), 1)
        columns = tsv_lines[0].split('\t')
        self.assertEqual(len(columns), 6)
        # Columns: normalized_id, '0', '0', '.', '.', 'Original name was ...'
        self.assertEqual(columns[1], '0')
        self.assertEqual(columns[2], '0')
        self.assertEqual(columns[3], '.')
        self.assertEqual(columns[4], '.')


if __name__ == '__main__':
    unittest.main()
