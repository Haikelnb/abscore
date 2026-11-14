"""Real antibody sequences used as test fixtures.

Every sequence here is copied from a public, MIT-licensed source -- never typed
from memory. Provenance for each is recorded below so a reviewer can check it.

Source
------
Merck/BioPhi-2021-publication, `data/tasks/humanness/imgt_mab_db.fa` and
`imgt_mab_db.tsv` (MIT License, Copyright 2021 Merck Sharp & Dohme Corp.),
which were in turn derived from IMGT/mAb-DB.

  https://github.com/Merck/BioPhi-2021-publication
  Prihoda D et al. BioPhi. MAbs 2022;14(1):2020203.

Retrieved 2026-11-12. FASTA records are keyed by INN number; `_seq1` and `_seq2`
are the two chains of a pair, in the order the file lists them (which is not
always heavy first -- see ADALIMUMAB below).
"""

# --- trastuzumab (INN 7637), labelled "Humanized" in imgt_mab_db.tsv ----------
# 7637_seq1
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
# 7637_seq2 -- kappa light chain
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

# --- rituximab (INN 7609), labelled "Chimeric" --------------------------------
# A chimeric antibody's variable domains are entirely murine; only the constant
# regions are human. So on V-domain humanness it is expected to score with the
# murine antibodies, leaving murine and humanized indistinguishable. Asserted.
RITUXIMAB_VH = (
    "QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAIYPGNGDTSYNQKFKG"
    "KATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWGAGTTVTVSA"
)
RITUXIMAB_VL = (
    "QIVLSQSPAILSASPGEKVTMTCRASSSVSYIHWFQQKPGSSPKPWIYATSNLASGVPVRFSGSGS"
    "GTSYSLTISRVEAEDAATYYCQQWTSNPPTFGGGTKLEIK"
)

# --- adalimumab (INN 7860), labelled "Human" ----------------------------------
# NOTE the file lists the light chain first for this entry. Kept in source order
# with explicit names so the ordering assumption is never silent.
ADALIMUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSG"
    "SGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK"
)
ADALIMUMAB_VH = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEG"
    "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)

# --- pembrolizumab (INN 9798), labelled "Humanized" ---------------------------
# Worth keeping as a fixture: ANARCI's top numbering HMM for this VH is
# *mouse* (bitscore 184.1) rather than human (180.5), despite the antibody being
# a marketed humanized therapeutic. The ANARCI species call is a numbering
# nearest-neighbour; this sequence shows why humanness is scored separately.
PEMBROLIZUMAB_VL = (
    "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARF"
    "SGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK"
)
PEMBROLIZUMAB_VH = (
    "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKN"
    "RVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
)

# --- abagovomab (INN 8659), labelled "Mouse" ----------------------------------
ABAGOVOMAB_VH = (
    "QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKG"
    "KATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS"
)

# --- elipovimab (INN 10988), labelled "Human" ---------------------------------
# A 26-residue CDR3, the longest in the panel, carrying 13 IMGT insertion codes:
#   105 ... 111 111A 111B 111C 111D 111E 111F 112G 112F 112E 112D 112C 112B
#   112A 112 113 ... 117
# Note that insertions at 112 run DESCENDING and sit BEFORE position 112. This
# is the IMGT "outwards from the top of the loop" convention, and it is the
# regression fixture for imgt_sort_key: a naive (number, letter) sort transposes
# these residues silently.
ELIPOVIMAB_VH = (
    "QMQLQESGPGLVKPSETLSLTCSVSGASISDSYWSWIRRSPGKGLEWIGYVHKSGDTNYNPSLKSR"
    "VHLSLDTSKNQVSLSLTGVTAADSGKYYCARTLHGRRIYGIVAFNEWFTYFYMDVWGTGTQVTVSS"
)
ELIPOVIMAB_VH_CDR3 = "ARTLHGRRIYGIVAFNEWFTYFYMDV"

# --- a non-antibody protein, for negative tests -------------------------------
# Human serum albumin (UniProt P02768) mature N-terminal region, residues 1-65.
# Source: UniProt P02768, retrieved 2026-11-12.
ALBUMIN_FRAGMENT = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS"

# A V-domain followed by a non-antibody tail. Numbering must trim to the
# V-domain, exactly as it must for a full-length chain carrying CH1-CH3.
VH_WITH_TAIL = TRASTUZUMAB_VH + ALBUMIN_FRAGMENT

# Origin labels as given in imgt_mab_db.tsv, for tests that assert ordering.
ORIGIN_LABELS = {
    "trastuzumab": "Humanized",
    "rituximab": "Chimeric",
    "adalimumab": "Human",
    "pembrolizumab": "Humanized",
    "abagovomab": "Mouse",
}

# --- abelacimab (INN 10881) lambda light chain, labelled "Human" --------------
# A lambda chain, closing the gap that the earlier fixtures left: they were all
# heavy or kappa. Confirms that lambda, like kappa, ends at IMGT 127 and carries
# J-PHE (F) rather than J-TRP at position 118. 10881_seq1.
ABELACIMAB_VL_LAMBDA = (
    "QSVLTQPPSASGTPGQRVTISCSGSSSNIGSNDVSWYQQLPGTAPKLLIYKNYNRPSGVPDRFSGSK"
    "SGTSASLAISGLQSEDEADYYCSAWDQRQFDVVFGGGTKLTVL"
)
