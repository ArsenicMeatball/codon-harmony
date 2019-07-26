#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony` package."""


import unittest

from codon_harmony import codon_harmony
from codon_harmony.codon_harmony import DictAsArg


class TestCodon_tools(unittest.TestCase):
    """Tests for `codon_tools` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.args_to_parse = [
            "--input",
            "misc/INPUT_LIST.fasta",
            "--output",
            "out.fasta",
        ]

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_args(self):
        """Test args parsing."""
        parser = codon_harmony.get_parser()
        parsed_args = parser.parse_args(self.args_to_parse)
        assert parsed_args.input == "misc/INPUT_LIST.fasta"
        assert parsed_args.output == "out.fasta"
        assert parsed_args.cycles == 10
        assert parsed_args.host == "413997"
        assert parsed_args.host_threshold == 0.1
        assert parsed_args.local_host_profile is None
        assert parsed_args.inner_cycles == 10
        assert parsed_args.local_homopolymer_threshold == 4
        assert parsed_args.max_relax == 0.1
        assert parsed_args.restriction_enzymes == [
            "NdeI",
            "XhoI",
            "HpaI",
            "PstI",
            "EcoRV",
            "NcoI",
            "BamHI",
        ]
        assert parsed_args.verbose == 0
        assert parsed_args.splice_sites
        assert parsed_args.start_sites

    def test_args_dict(self):
        """Test dictionary passed as an argument"""

        """Test for default arguments"""
        test_dict = {"input": "input.fasta"}
        parsed_args = DictAsArg(test_dict)
        assert parsed_args.input == "input.fasta"
        assert parsed_args.output == "out.fasta"
        assert parsed_args.cycles == 10
        assert parsed_args.host == "413997"
        assert parsed_args.host_threshold == 0.1
        assert parsed_args.local_host_profile is None
        assert parsed_args.inner_cycles == 10
        assert parsed_args.local_homopolymer_threshold == 4
        assert parsed_args.max_relax == 0.1
        assert parsed_args.restriction_enzymes == [
            "NdeI",
            "XhoI",
            "HpaI",
            "PstI",
            "EcoRV",
            "NcoI",
            "BamHI"
        ]
        assert parsed_args.verbose == 0
        assert not parsed_args.splice_sites
        assert not parsed_args.start_sites
        assert not parsed_args.one_line_fasta
        assert parsed_args.run

        """Test for full dict"""
        test_dict = {
            "input": "misc/INPUT_LIST.fasta",
            "host": "413998",
            "host_threshold": 0.2,
            "local_homopolymer_threshold": 5,
            "cycles": 12,
            "inner_cycles": 12,
            "max_relax": 0.2,
            "restriction_enzymes": "NdeI XhoI HpaI PstI EcoRV NcoI".split(),
            "splice_sites": True,
            "start_sites": True,
            "verbose": 1,
            "one_line_fasta": True,
            "output": "output.fasta",
            "run": False,
        }
        parsed_args = DictAsArg(test_dict)
        assert parsed_args.input == "misc/INPUT_LIST.fasta"
        assert parsed_args.output == "output.fasta"
        assert parsed_args.cycles == 12
        assert parsed_args.host == "413998"
        assert parsed_args.host_threshold == 0.2
        assert parsed_args.local_host_profile is None
        assert parsed_args.inner_cycles == 12
        assert parsed_args.local_homopolymer_threshold == 5
        assert parsed_args.max_relax == 0.2
        assert parsed_args.restriction_enzymes == [
            "NdeI",
            "XhoI",
            "HpaI",
            "PstI",
            "EcoRV",
            "NcoI"
        ]
        assert parsed_args.verbose == 1
        assert parsed_args.splice_sites
        assert parsed_args.start_sites
        assert parsed_args.one_line_fasta
        assert not parsed_args.run

        """Integrity test - dict passing through the code"""
        codon_harmony.main({"input": "misc/INPUT_LIST.fasta"})

    def test_main(self):
        """Test `codon_harmony` main function"""
        # all this does is make sure that nothing crashes --
        # more of an integration test than a unit test with
        # each unit being separately tested.

        # test three sequences:
        #   one that can be optimized
        #   one that cannot be optimized (given max_relax = 0.1)
        #   one that triggers a warning about GC content
        codon_harmony.main(self.args_to_parse)
        assert True

    def test_harmonize_sequence(self):
        """Test `codon_harmony._harmonize_sequence` function"""
        # make sure the DNA is correct and labels are preserved
        from Bio.Alphabet import IUPAC
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from codon_harmony.util import codon_use
        from codon_harmony.data import RestrictionEnzymes

        parser = codon_harmony.get_parser()
        args = parser.parse_args(self.args_to_parse)

        cut, hp, cra = codon_use.host_codon_usage(args.host, args.host_threshold)
        rest_enz = RestrictionEnzymes(args.restriction_enzymes)

        seq_records = []
        seq_records.append(
            SeqRecord(
                Seq("HHHHHHHHHH", IUPAC.protein),
                id="test_sequence1",
                name="Test1",
                description="can be optimized with `max_relax` set to 0.1",
            )
        )

        seq_records.append(
            SeqRecord(
                Seq("ACDEFGHIKLMNPQRSTVWY", IUPAC.protein),
                id="test_sequence2",
                name="Test2",
                description="cannot be optimized with `max_relax` set to 0.1",
            )
        )

        seq_records.append(
            SeqRecord(
                Seq("FFFFFFFFFFFF", IUPAC.protein),
                id="test_sequence3",
                name="Test3",
                description="can be optimized with `max_relax` set to 0.1, has extreme GC content",
            )
        )

        for seq_record in seq_records:
            dna_sequence = codon_harmony._harmonize_sequence(
                seq_record, args, cut, hp, cra, rest_enz
            )

        label, seq = dna_sequence["dna"].strip().split("\n")
        assert label == ">{} {}".format(seq_record.id, seq_record.description)
        assert Seq(seq, IUPAC.unambiguous_dna).translate() == seq_record.seq
