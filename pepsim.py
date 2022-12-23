#!/usr/bin/env python
# AUTHOR: WSN
# CREATE DATE: 4 May 2017
import sys
import peptides

USAGE = """USAGE: computePeptideSimilarities.py <file1> <file2>

  Compute pairwise similarities among two sets of peptides.  Each
  input file contains one peptide sequence per line, with
  modifications indicated by brackets.  Output is a three-column,
  tab-delimited file with peptide1, peptide2, and score. The score is
  simply the proportion of matched b- and y-ions.  If an m/z threshold
  is specified, then four additional columns (mz1, mz2, mz
  difference, ppm difference) are printed.

  Options:

    --mz-thresh <value>   Threshold (in ppm) for computing similarities.
                          With this option, the input must have a second
                          column containing the peptide mass.
                          Default = no threshold.

    --frag-bin-size <value> Bin size used on fragment m/z axis.
                            Default = 0.05.

    --min-score <value>    Minimum score needed for which result is
                           store. Default = 0.25.

    --static-mods <string> A comma-separated list of static mods,
                           formatted like <amino acid>:<mass shift>.
                           Use "nterm" or "cterm" to indicate terminal
                           modifications.  C+57.02146 is added to
                           every set of mods, unless some other C
                           modification is included.

    --skip-same    Don't bother to compute the similarity between two
                   identical peptides.

"""

###############################################################################
def parse_static_mods(my_string):
  """
  Parse a static mod string (see USAGE) into a dictinoary.
  Key = amino acid, value = mass offset
  """
  return_value = {}
  
  for one_mod in my_string.split(","):
    words = one_mod.split(":")
    return_value[words[0]] = float(words[1])

  # Add in carbamidomethylation.
  if ("C" not in return_value):
    return_value["C"] = 57.02146
  
  return(return_value)
  
  
###############################################################################
def make_binidx_matchcount_map(mzs, fragment_min_mz, fragment_max_mz, bin_size):
    """
    Utility function for calc_proportion_fragments_incommon.
    Notionally bin the mzs from a list of mzs using the specified bin
    size (don't actually build the binned array), with specified lower
    and upper limits. Return a map from bin indexes to a count of
    fragments
    :param mzs: 
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:  a map from bin index to a count of fragments in that bin
    """
    binidx_matchcount_map = {}
    for mz in mzs:
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx not in binidx_matchcount_map:
            binidx_matchcount_map[bin_idx] = 0
        binidx_matchcount_map[bin_idx] += 1
    return binidx_matchcount_map

###############################################################################
def calc_proportion_fragments_incommon(pepseq1, peptide1_mods, nterm1, cterm1,
                                       pepseq2, peptide2_mods, nterm2, cterm2,
                                       binsize):
    """
    Determine all the fragment mzs for each peptide. Bin the fragments from
    each peptide.  Calculate the fraction of fragments that fall into a bin
    with a fragment from the other peptide.

    :param pepseq1: First peptide sequence
    :param peptide1_mods: List of mass offsets, same length as peptide.
    :param nterm1: N-terminal modification mass of first peptide.
    :param cterm1: C-terminal modification mass of first peptide.
    :param pepseq2: Second peptide sequence
    :param peptide2_mods: List of mass offsets, same length as peptide.
    :param nterm2: N-terminal modification mass of second peptide.
    :param cterm2: C-terminal modification mass of second peptide.
    :param binsize: Size of m/z bins.
    :return: Fraction of matched peaks.
    """

    # Compute fragments and put them in bins.
    mzs_1 = peptides.calc_theoretical_peak_mzs(pepseq1, [1], peptide1_mods, 
                                               200, 3000, nterm1, cterm1)
    mzs_2 = peptides.calc_theoretical_peak_mzs(pepseq2, [1], peptide2_mods,
                                               200, 3000, nterm2, cterm2)
    bin_count_map_1 = make_binidx_matchcount_map(mzs_1, 200, 3000, binsize)
    bin_count_map_2 = make_binidx_matchcount_map(mzs_2, 200, 3000, binsize)

    # Count matched bins.
    n_fragments_in_matched_bins = 0
    for binidx in bin_count_map_1:
        if binidx in bin_count_map_2:
            n_fragments_in_matched_bins += (bin_count_map_1[binidx]
                                            + bin_count_map_2[binidx])

    return float(n_fragments_in_matched_bins) / (len(mzs_1) + len(mzs_2))


###############################################################################
def sort_peptide_list(peptide_list, mass_list):
  # Create a list of tuples.
  tuples = list(zip(peptide_list, mass_list))
  # Sort the list by the second element in each.
  tuples.sort(key=lambda tup: tup[1])
  # Undo the zip operation.
  return(zip(*tuples))

###############################################################################
def read_peptides(in_file_name):
  """
  Read peptides from the first tab-delimited column of a given file.
  """
  peptide_list = []
  mass_list = []

  inFile = open(in_file_name, "r")
  # for new tide-index (combined target and decoy peptie files)
  header = inFile.readline()
  for line in inFile:
    words = line.rstrip().split("\t")
    peptide_list.append(words[0])
    if (len(words) > 1):
      mass_list.append(float(words[2]))
  inFile.close()
  if (len(peptide_list) != len(mass_list)):
    sys.stderr.write("Error: Found %d peptides but only %d masses.\n"
                     % (len(peptide_list), len(mass_list)))
    sys.exit(1)
  sys.stderr.write("Read %d peptides from %s.\n"
                   % (len(peptide_list), in_file_name))

  # Sort the peptides by mass.
  (peptide_list, mass_list) = sort_peptide_list(peptide_list, mass_list)

  return(peptide_list, mass_list)


###############################################################################
def parse_mods(pepseq_with_mods, static_mods):
    """
    Parse a modified peptide sequence string.

    :param pepseq_with_mods: Peptide string with interpolated bracketed mods.  
    :param static_mods: Dictionary of static mods. 
                        Key = amino acid, value = mass offset.
    :return: A list of amino acids, a list of modification values, and 
             the n-term and c-term terminal modification values.
    """
    aa_seq_list = []
    modifications = []
    nterm_delta = 0.0
    cterm_delta = 0.0

    aaseq_position = -1
    modpepseq_position = 0
    while modpepseq_position < len(pepseq_with_mods):
        my_char = pepseq_with_mods[modpepseq_position]
        if my_char.isalpha():

            # Create an unmodified amino acid and add it to the growing list.
            aa_seq_list.append(my_char)
            modifications.append(0.0)
            aaseq_position += 1
            modpepseq_position += 1
        elif ( my_char == '[' ):
            end_pos = (pepseq_with_mods[modpepseq_position + 1:].index(']') 
                       + modpepseq_position + 1)
            mod_mass = float(pepseq_with_mods[modpepseq_position + 1:end_pos])

            # Store a modification associated with the current position.
            modifications[aaseq_position] = mod_mass
            modpepseq_position = end_pos + 1
        else:
            sys.stderr.write("Invalid character (%s) at position %d.\n"
                             % (my_char, modpepseq_position))
            sys.exit(1)

    # Add in the static mods.
    for index in range(0, len(aa_seq_list)):
        amino = aa_seq_list[index]
        if (amino in static_mods):
          modifications[index] += static_mods[amino]
    if ("nterm" in static_mods):
        nterm_delta = static_mods["nterm"]
    if ("cterm" in static_mods):
        cterm_delta = static_mods["cterm"]

    return(aa_seq_list, modifications, nterm_delta, cterm_delta)

###############################################################################
def ppm_difference (mass1, mass2):
  """
  Calculate the difference in parts-per-million between two peptide masses.
  :param mass1: Mass of the first peptide.
  :param mass2: Mass of the second peptide.
  :return: Difference in ppm.
  """

  return(1e6 * (mass2 - mass1) / (0.5 * (mass1 + mass2)))
  
###############################################################################
def main():
  global USAGE

  # Set default values for parameters.
  frag_bin_size = 0.05 # Size of fragment m/z bins in Da.
  mz_thresh = 0.0 # Threshold in ppm below which similarities are computed.
  min_score = 0.25 # Score threshold for which result is kept
  static_mods = {"C": 57.02146} # Key = amino acid, value = mass offset
  skip_same = False

  # Parse the command line.
  sys.argv = sys.argv[1:]
  while (len(sys.argv) > 2):
    next_arg = sys.argv[0]
    sys.argv = sys.argv[1:]
    if (next_arg == "--frag-bin-size"):
      frag_bin_size = float(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "--mz-thresh"):
      mz_thresh = float(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "--min-score"):
      min_score = float(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "--static-mods"):
      static_mods = parse_static_mods(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "--skip-same"):
      skip_same = True
    else:
      sys.stderr.write("Invalid option (%s).\n" % next_arg)
      sys.exit(1)
  if (len(sys.argv) != 2):
    sys.stderr.write(USAGE)
    sys.exit(1)
  in_file_name1 = sys.argv[0]
  in_file_name2 = sys.argv[1]

  # Read in all the peptides.
  (list1, masses1) = read_peptides(in_file_name1)
  (list2, masses2) = read_peptides(in_file_name2)

  start_index = 0

  # Compare each pair of peptides.
  num_pairs = 0
  num_printed = 0
  for index1 in range(0, len(list1)):
    peptide1 = list1[index1]
    (unmodified_peptide1, peptide1_mods, nterm1, cterm1) \
       = parse_mods(peptide1, static_mods)

    if ((index1+1) % 1000 == 0):
      print(f"{(index1 / len(list1)) * 100:.2f}% complete. " +
            f"Band width = {index2 - start_index} " +
            f"({((index2 - start_index) / len(list2)) * 100:.2f}%).",
            file=sys.stderr)

    for index2 in range(start_index, len(list2)):
      peptide2 = list2[index2]
      num_pairs += 1

      mass_diff = ppm_difference(masses1[index1], masses2[index2])

      # If the second mass is outside the range, terminate this loop.
      if (mz_thresh != 0.0) and (mass_diff > mz_thresh):
        break

      # If we haven't gotten to the range yet, increment the start index.
      if (mz_thresh != 0.0) and (mass_diff < -1.0 * mz_thresh):
        start_index = index2
        continue

      # Don't bother if they are the same peptide.
      if (skip_same and
          (peptide1.replace("I", "L") == peptide2.replace("I", "L"))):
        continue

      # Don't bother computing similarity if the masses are too different.
      if ((mz_thresh != 0.0) and (abs(mass_diff) <= mz_thresh)):

        (unmodified_peptide2, peptide2_mods, nterm2, cterm2) \
         = parse_mods(peptide2, static_mods)

        similarity = calc_proportion_fragments_incommon(
          unmodified_peptide1, peptide1_mods, nterm1, cterm1,
          unmodified_peptide2, peptide2_mods, nterm2, cterm2,
          frag_bin_size)

        if similarity >= min_score:
          num_printed += 1
          print(f"{peptide1}\t{peptide2}\t{similarity:.4g}", end="")
          if (mz_thresh != 0.0):
            print(f"\t{masses1[index1]:.4f}" +
                  f"\t{masses2[index2]:.4f}" +
                  f"\t{masses1[index1] - masses2[index2]:.4f}" +
                  f"\t{mass_diff:.4f}", end="")
          print("")

  print(f"Evaluated {num_pairs} out of {len(list1) * len(list2)} possible " +
        f"pairs ({(num_pairs / (len(list1)*len(list2))) * 100:.2f})",
        file=sys.stderr)
  print(f"Printed {num_printed} out of {num_pairs} evaluated pairs " +
        f"({(num_printed / num_pairs) * 100:.2f}).", file=sys.stderr)


if __name__ == "__main__":
  main()
