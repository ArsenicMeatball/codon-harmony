from dataclasses import dataclass, field


@dataclass
class DictAsArg:
    """
    allows users to pass a dictionary as an argument
    example of usage
    the_dictionary = {'host': 'my_host_id', 'verbose': 2}
    codon_harmony.runner(the_dictionary)
    """
    input: str = "input.fasta"
    host: str = "413997"
    host_threshold: float = 0.1
    cycles: int = 10
    local_homopolymer_threshold: int = 4
    inner_cycles: int = 10
    max_relax: float = 0.1
    longline = "NdeI XhoI HpaI PstI EcoRV NcoI BamHI"
    restriction_enzymes: [str] = field(default_factory=longline.split)
    splice_sites: bool = False
    start_sites: bool = False
    local_host_profile: str = None
    verbose: int = 0
    one_line_fasta: bool = False
    output: str = "out.fasta"
    run: bool = True

    @classmethod
    def from_dict(cls, dictionary):
        inst = cls()
        inst.__dict__.update(
            {k: v for k, v in dictionary.items() if k in inst.__dict__}
        )
        return inst
