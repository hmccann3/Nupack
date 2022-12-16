from nupack import *
import pandas as pd
from string import ascii_lowercase
import Levenshtein as lev
import re
### Define physical model for NUPACK. Make sure to adjust to physiological temperature
my_model = Model(material="rna", celsius=28)
### Input file containing anchors and compactors
orig_df = pd.read_csv("compactor_summary.tsv", sep="\t")
struct_lev_dict = {}
j = 0
### Only include compactors above average compactor length
compactor_length = int(orig_df["compactor_valid"].apply(len).mean())
orig_df = orig_df[orig_df["compactor_valid"].str.len() >= compactor_length]
###Iterate through anchors
for anchor in orig_df.anchor.unique():
    j += 1
    ### Only run nupack when there are more than five compactors for an anchor
    if len(orig_df[orig_df["anchor"] == anchor]) > 5:
        df = orig_df[orig_df["anchor"] == anchor]
        alphabet = []
        for c in ascii_lowercase:
            for c_1 in ascii_lowercase:
                alphabet.append(c + "_" + c_1)
        i = 0
        compactor_dict = {}
        a_dict = {}
        swap_dict = {"A": "U", "T": "A", "G": "C", "C": "G"}
        for compactor in df["compactor_valid"]:
            ### Standardize compactor length and convert from cDNA to RNA
            short_compactor = compactor[:compactor_length]
            rna_compactor = ""
            for char in short_compactor:
                if char in swap_dict:
                    rna_compactor += swap_dict[char]
                else:
                    print("Skipping base N")
            ###Calculate compactor weights based on valid_local_proportion from compactor_summary file 
            valid_proportion = df[df["compactor_valid"] == compactor]["valid_local_proportion"].astype("float").iloc[0]
            compactor_weighting = valid_proportion/(df["valid_local_proportion"].sum())
            compactor_dict["(" + alphabet[i] + ")"] = [rna_compactor, compactor, compactor_weighting]
            a_dict[Strand(rna_compactor, name=alphabet[i])] = 5e-6
            i += 1
        ###Run nupack on all compactors for this anchor 
        t1 = Tube(strands=a_dict, complexes=SetSpec(max_size=1), name="Tube t1")
        tube_result = tube_analysis(
            tubes=[t1],
            model=my_model,
            compute=["mfe", "subopt"],
            options={"energy_gap": 1.5},
        )
        ###For complex analysis
        """
        for result in tube_result:
            walker_result = tube_result[result]
            for i, s in enumerate(walker_result.mfe):
                f.write('>'+result.name + '\n')
                f.write(compactor_dict[result.name] + '\n')
                f.write('%s'%(s.structure) + '\n')
                #print('    %2d: %s (%.2f kcal/mol)' % (i, s.structure, s.energy))
        f.close()
        """
        ###For tube analysis
        struct_dict = {}
        for key in compactor_dict:
            walker_result = tube_result[key]
            subopt = walker_result.subopt
            for i, s in enumerate(walker_result.mfe):
                structs = [item[0] for item in subopt]
                # energies = [item[1] for item in subopt]
                ###Get second best structure from suboptimal structures
                unique_structs = [*set(structs)]
                if len(unique_structs) > 1:
                    second_best_struct = unique_structs[1]
                else:
                    second_best_struct = unique_structs[0]
                struct_dict[compactor_dict[key][0]] = [
                    "%s" % (s.structure),
                    s.energy,
                    compactor_dict[key][0],
                    compactor_dict[key][1],
                    "%s" % (second_best_struct),
                    compactor_dict[key][2]
                ]
        known_comparisons = {}
        ###Calculate levenshtein distances and covariance for compactors/structures
        for key, val in struct_dict.items():
            struct = val[0]
            mfe = val[1]
            short_compactor = val[2]
            full_compactor = val[3]
            second_best_struct = val[4]
            compactor_weight = val[5]
            struct_lev_dist_total = 0
            subopt_struct_lev_dist_total = 0
            lev_dist_total = 0
            covarying_pairs_total = 0
            for val2 in struct_dict.values():
                compactor2_weight = val2[5]
                weight = compactor_weight * compactor2_weight
                if val[3] != val2[3]:
                    if (val2[3], val[3]) in known_comparisons:
                        entry = known_comparisons[(val2[3], val[3])]
                        struct_lev_dist = entry[0]
                        subopt_struct_lev_dist = entry[1]
                        seq_lev_dist = entry[2]
                        covarying_pairs = entry[3]
                        known_comparisons.pop(val2[3], val[3])
                    else:
                        struct2 = val2[0]
                        second_best_struct2 = val2[4]
                        short_compactor2 = val2[2]
                        ###Get lowest levenshtein distance between top two structures of two compactors
                        struct_lev_dist = lev.distance(struct, struct2)
                        subopt_struct_lev_dist = min(
                            [
                                struct_lev_dist,
                                lev.distance(struct, second_best_struct2),
                                lev.distance(second_best_struct, struct2),
                                lev.distance(second_best_struct, second_best_struct2),
                            ]
                        )
                        seq_lev_dist = lev.distance(short_compactor, short_compactor2)
                        mismatches = [i for i in range(len(short_compactor)) if short_compactor[i] != short_compactor2[i]]
                        ###Count number of covarying pairs between two compactors
                        covarying_pairs = 0
                        for pos in mismatches:
                            if struct[pos] == struct2[pos] and struct[pos] == '(':
                                left_count = 1
                                right_count = 0
                                i = pos + 1
                                while i < len(struct):
                                    if struct[i] == '(':
                                        left_count += 1
                                    elif struct[i] == ')':
                                        right_count += 1
                                    if left_count == right_count:
                                        if i in mismatches:
                                            covarying_pairs += 1
                                        break
                                    i += 1
                        known_comparisons[(val[3], val2[3])] = [struct_lev_dist, subopt_struct_lev_dist, seq_lev_dist, covarying_pairs]
                    subopt_struct_lev_dist_total += subopt_struct_lev_dist * weight
                    struct_lev_dist_total += struct_lev_dist * weight
                    # struct_lev_dist_total += struc_distance(struct, struct2)
                    lev_dist_total += seq_lev_dist * weight
                    # lev_dist_total += seq_distance(short_compactor, short_compactor2)
                    covarying_pairs_total += covarying_pairs * weight
            struct_lev_dict[full_compactor] = [
                struct_lev_dist_total / (len(struct_dict) - 1),
                lev_dist_total / (len(struct_dict) - 1),
                struct,
                mfe,
                mfe * compactor_weight,
                subopt_struct_lev_dist_total / (len(struct_dict) - 1),
                covarying_pairs_total / (len(struct_dict) - 1)
            ]
newdf = pd.DataFrame.from_dict(
    struct_lev_dict,
    orient="index",
    columns=[
        "avg_pairwise_lev_struct_dist",
        "avg_pairwise_lev_dist",
        "nupack_structure",
        "mfe",
        "weighted_mfe",
        "avg_pairwise_best_subopt_lev_struct_dist",
        "avg_pairwise_covarying_pairs"
    ],
)
newdf.reset_index(inplace=True)
newdf = newdf.rename(columns={"index": "compactor_valid"})
final_df = orig_df.merge(newdf, how="left", on="compactor_valid")
final_df.dropna(subset=["avg_pairwise_lev_struct_dist"], inplace=True)
final_df["lev_dist_struct_dist_ratio"] = (
    final_df["avg_pairwise_lev_struct_dist"] / final_df["avg_pairwise_lev_dist"]
)
final_df["best_subopt_lev_dist_struct_dist_ratio"] = (
    final_df["avg_pairwise_best_subopt_lev_struct_dist"] / final_df["avg_pairwise_lev_dist"]
)
final_df["weighted_mfe_mean"] = (final_df["weighted_mfe"].groupby(final_df["anchor"]).transform("sum"))
final_df["lev_dist_struct_dist_mean"] = (
    final_df["lev_dist_struct_dist_ratio"].groupby(final_df["anchor"]).transform("mean")
)
final_df["best_subopt_lev_dist_struct_dist_mean"] = (
    final_df["avg_pairwise_best_subopt_lev_dist_struct_dist_ratio"]
    .groupby(final_df["anchor"])
    .transform("mean")
)
final_df["covarying_pairs_mean"] = (
    final_df["avg_pairwise_covarying_pairs"]
    .groupby(final_df["anchor"])
    .transform("mean")
)
###Sort output
# final_df.sort_values('mfe_mean', ascending = True, inplace = True)
final_df.sort_values(
    "best_subopt_lev_dist_struct_dist_mean", ascending=True, inplace=True
)
final_df = final_df[final_df.columns.drop(list(final_df.filter(regex="Unnamed")))]
final_df.drop("weighted_mfe", inplace=True, axis=1)
###Output file with anchors, compactors and nupack structures and scores
final_df.to_csv("compactor_summary_nupack.tsv", sep="\t")
