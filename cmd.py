from chimerax.core.commands import CmdDesc
from chimerax.atomic import AtomsArg
# imports for sequence alignment
from chimerax.alignment_algs import SmithWaterman
from chimerax import sim_matrices
import csv
hello_world_desc = CmdDesc()


class DistanceCalculator:
    def __init__(self, session):
        self.session = session
        self.distances = []

# Read in FASTA file

        # Read in fasta file
        def read_fasta_file(fasta_file):
            """
            :param fasta_file: str Path to the input FASTA file.
            :return:  dict
                 Dictionary where keys are UniProt IDs and values are the corresponding amino acid sequences.
            """

            uniprot_to_seq = {}
            current_id = None
            current_seq = ""

            with open(fasta_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_id and current_seq:
                            uniprot_to_seq[current_id] = current_seq
                        header_parts = line[1:].split("|")
                        if len(header_parts) >= 2:
                            current_id = header_parts[1]
                        else:
                            current_id = line[1:]
                        current_seq = ""
                    else:
                        current_seq += line

            if current_id and current_seq:
                uniprot_to_seq[current_id] = current_seq

            return uniprot_to_seq
# Read in XL-MS dataset
        def read_xlms_data(tsv_file):
            """
            :param tsv_file: dataset containing results from XL-MS data
            :return: dict containing positions of crosslink pairs
            """


            accession_to_pos = []

            with open(tsv_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    acc_a = row['Accession A']
                    pos_a = int(row['Position A'])
                    acc_b = row['Accession B']
                    pos_b = int(row['Position B'])
                    accession_to_pos.append((acc_a, pos_a, acc_b, pos_b))

            return accession_to_pos

        data = read_xlms_data("/Users/catherinenguyen/Downloads/media-2/evidence1.tsv")
        pdb = self.session.models
        file = read_fasta_file("/Users/catherinenguyen/Downloads/P02671.fasta")

        def find_best_alignment(uniprot_to_seq, models):
            """

            :param dict uniprot_to_seq: Dictionary of FASTA file including UniProt accessions and full sequences
            :param models: structural models chain sequences on the ChimeraX PDB which is to be aligned with Uniprot sequences
            :return dict: best matched sequences including:
                        Uniprot sequence, uniprot aligned sequence, pdb aligned sequence, model id,
                        chain id, pdb start position, alignment score.
            """
            best_matches = {}

            for protein_id, sequence in uniprot_to_seq.items():
                best_score = -1
                best_result = None

                for model in models:
                    model_id = model.id_string
                    if hasattr(model, 'chains') and model.chains:
                        for chain in model.chains:
                            chain_seq = chain.characters
                            pdb_start_pos = chain.numbering_start

                            # Use Smith-Waterman alignment
                            matrix = sim_matrices.matrix("BLOSUM-62")
                            score, alignment = SmithWaterman.align(sequence, chain_seq, matrix, 11, 1)
                            uniprot_aligned, pdb_aligned = alignment

                            # Keep only the best alignment
                            if score > best_score:
                                best_score = score
                                best_result = {
                                    "sequence": sequence,
                                    "aligned_ref": uniprot_aligned,
                                    "aligned_chain": pdb_aligned,
                                    "model": model_id,
                                    "chain": chain.chain_id,
                                    "pdb_start": pdb_start_pos,
                                    "score": score
                                }

                # Save best result for this protein
                if best_result:
                    best_matches[protein_id] = best_result

            return best_matches

        results = find_best_alignment(file, pdb)
        print(str(results))

        def map_crosslinks_to_pdb(tsv_data, alignments):
            """

            :param tsv_data: dict containing positions of crosslink pairs
            :param alignments: dict containing best matched sequences
            :return: mapped uniprot residue positions to PDB positions in ChimeraX
            """

            mapped_coords = []

            for acc_a, pos_a, acc_b, pos_b in tsv_data:
                coords = []
                for acc, pos in [(acc_a, pos_a), (acc_b, pos_b)]:
                    if acc not in alignments:
                        coords.append(None)
                        continue

                    aln = alignments[acc]
                    chain_aligned = aln["aligned_chain"]
                    seq_uniprot = aln["aligned_ref"]

                    if len(chain_aligned) < 6:
                        coords.append(None)
                        continue

                    # Find where aligned region starts in UniProt
                    aligned_region = chain_aligned.replace("-", "")
                    uniprot_start = seq_uniprot.find(aligned_region) + 1

                    # Map UniProt residue pos to PDB pos
                    pdb_start_pos = aln["pdb_start"]
                    difference = pos - uniprot_start
                    pdb_res_pos = difference + pdb_start_pos

                    if difference < len(aligned_region):
                        coord = "#" + str(aln['model']) + "/" + aln['chain'] + ":" + str(pdb_res_pos) + "@ca"

                    else:
                        coord = None

                    coords.append(coord)
                # Store as tuple (coord_a, coord_b)
                mapped_coords.append(coords)

            return mapped_coords

        global mapped_global
        mapped_global = map_crosslinks_to_pdb(data, results)

    def calc_distances(self, a1, a2):
        """
        Describe me
        :param str a1: the first atom from which the crosslink is made
        :param str a2: the second atom from which the crosslink is made
        :return: the distance of the crosslink between the two atoms;
                visual pseudo bonds between the atom pairs with distance
        """
        # Make a new group
        pbg = self.session.pb_manager.get_group("crosslinks")

        if pbg.id is None:  # If the group is new, add it to the session
            self.session.models.add([pbg])

        pbg.display = True

        # Assign atom pairs to each atom
        a1 = AtomsArg.parse(a1, self.session)
        a2 = AtomsArg.parse(a2, self.session)
        pb = pbg.new_pseudobond(a1[0][0], a2[0][0])
        pb.radius = 0.1
        self.session.logger.info(str(pb.length))
        self.session.logger.info(str(pb))

        # Measure distance - shown on PDB
        distance = self.session.pb_dist_monitor.add_group(pbg)
        # Print distance
        self.session.logger.info.pb_dist_monitor.add_group(pbg)

def hello_world(session):   # Cant take mapped as a parameter so have to store mapped globally
    """

    :param session: ChimeraX open session
    :return: Distance measured between two atom pairs or warnings of
            missing structures or missing model.
    """
    calculator = DistanceCalculator(session)

    atom_pairs = []
    for a1, a2 in mapped_global:
        if a1 is not None and a2 is not None:
            atom_pairs.append((a1, a2))
        else:
            print(f"Skipping pair ({a1}, {a2}) due to error: Missing model")

    # print(atom_pairs)
    # Calculate distances for each atom pair
    for atom1, atom2 in atom_pairs:
        try:
            calculator.calc_distances(atom1, atom2)

        except Exception:
            print(f"Skipping pair ({atom1}, {atom2}) due to error: Missing structure")
            continue
