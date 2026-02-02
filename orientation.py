import argparse
import logging
import mappy as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

START_SEQ = (
    "AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATATCGCCAGCAGC"
    "ACATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGGCGGGCAGATACTTTAACCAA"
    "TATAGGAATACAAGACAGACAAATAAAAATGACAGAGTACACAACATCCATGAACCGCATCAGCACCACC"
    "ACCATTACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
)

def find_best_orientation_mappy(contig_seq):
    """Uses mappy (minimap2) to find the start sequence location."""
    # Index the contig (as a reference)
    # Using default presets (map-ont or map-hifi) is efficient for assembly mapping
    aligner = mp.Aligner(seq=str(contig_seq))
    if not aligner:
        return None, None

    # Map the START_SEQ (as a query) to the indexed contig
    hits = list(aligner.map(START_SEQ))
    
    if not hits:
        return None, None

    # Sort hits by mapping quality or score to get the best one
    best_hit = max(hits, key=lambda x: x.mlen)
    
    # mappy.Alignment.strand: +1 for forward, -1 for reverse
    orientation = "forward" if best_hit.strand == 1 else "reverse"
    
    # r_st is the start position on the reference (contig)
    return orientation, best_hit.r_st

def main():
    parser = argparse.ArgumentParser(description="Fast orientation using mappy.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-output", "--output", required=True, help="Output FASTA")
    args = parser.parse_args()

    oriented_records = []
    
    for record in SeqIO.parse(args.input, "fasta"):
        contig_len = len(record.seq)
        logger.info(f"Processing: {record.id} | Length: {contig_len} bp")
        
        orientation, start_pos = find_best_orientation_mappy(record.seq)
        
        if orientation == "forward":
            logger.info(f"  -> Match FORWARD at {start_pos}. Rotating.")
            rotated_seq = record.seq[start_pos:] + record.seq[:start_pos]
            oriented_records.append(SeqRecord(rotated_seq, id=record.id, description=record.description))
        elif orientation == "reverse":
            logger.info(f"  -> Match REVERSE at {start_pos}. Rotating.")
            rev_seq = record.seq.reverse_complement()
            # Calculate position on the RC strand
            # Note: best_hit.r_st is relative to the forward strand index
            # Rotating RC requires using the corresponding position
            rotated_seq = rev_seq[start_pos:] + rev_seq[:start_pos]
            oriented_records.append(SeqRecord(rotated_seq, id=record.id, description=record.description))
        else:
            logger.warning("  -> No match found. Keeping original.")
            oriented_records.append(record)

    SeqIO.write(oriented_records, args.output, "fasta")

if __name__ == "__main__":
    main()
