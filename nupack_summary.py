from nupack import *
import pandas as pd
from string import ascii_lowercase
import Levenshtein as lev

### Define physical model for nupack
#my_model = Model(material='rna', celsius=28)
my_model = Model(material='rna', celsius=37)
###Input file containing anchors and compactors
orig_df = pd.read_csv('compactor_summary.tsv', sep="\t")
struct_lev_dict = {}
lev_dict = {}
j = 0
for anchor in orig_df.anchor.unique()[:100]:
    j += 1
    #Only run nupack when there is more than one compactor
    if(len(orig_df[orig_df['anchor'] == anchor]) > 1):
        df= orig_df[orig_df['anchor'] == anchor]
        alphabet = []
        for c in ascii_lowercase:
            for c_1 in ascii_lowercase:
                alphabet.append(c + '_' + c_1)
        i = 0
        compactor_dict = {}
        a_dict = {}
        swap_dict = {'A': 'U', 'T':'A', 'G':'C', 'C':'G'}
        for compactor in df['compactor_valid']:
            rna_compactor = ''
            for char in compactor:
                if char in swap_dict:
                    rna_compactor += swap_dict[char]
                else:
                    print("Skipping base N")
            compactor_dict['('+alphabet[i]+')'] = compactor
            a_dict[Strand(rna_compactor, name=alphabet[i])] = 5e-6
            i+=1
        t1 = Tube(strands=a_dict, complexes=SetSpec(max_size=1), name='Tube t1')
        tube_result = tube_analysis(tubes=[t1], compute=['mfe'], model=my_model)
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
            for i, s in enumerate(walker_result.mfe):
                struct_dict[compactor_dict[key]] = ['%s'%(s.structure), s.energy, compactor_dict[key]]
        for key,val in struct_dict.items():
            struct = val[0]
            mfe = val[1]
            compactor = val[2]
            struct_lev_dist_total = 0
            lev_dist_total = 0
            for val2 in struct_dict.values():
                struct2 = val2[0]
                compactor2 = val2[2]
                struct_lev_dist_total += lev.distance(struct, struct2)
                lev_dist_total += lev.distance(compactor, compactor2)
            struct_lev_dict[key] = [struct_lev_dist_total/(len(struct_dict)), lev_dist_total/len(struct_dict), struct, mfe]
newdf = pd.DataFrame.from_dict(struct_lev_dict, orient='index', columns=['lev_struct_dist', 'avg_pairwise_lev_dist', 'nupack_structure', 'mfe'])
newdf.reset_index(inplace=True)
newdf = newdf.rename(columns = {'index':'compactor_valid'})
final_df = orig_df.merge(newdf, how='left', on='compactor_valid')
final_df.dropna(subset = ['lev_struct_dist'], inplace = True)
final_df['lev_dist_struct_dist_ratio'] = final_df['avg_pairwise_lev_dist']/final_df['lev_struct_dist']
final_df['mfe_mean'] = final_df['mfe'].groupby(final_df['anchor']).transform('mean')
###Sort by minimum free energy
final_df.sort_values('mfe_mean', ascending = True, inplace = True)
final_df = final_df[final_df.columns.drop(list(final_df.filter(regex='Unnamed')))]
###Output file with anchors, compactors and nupack structures and scores
final_df.to_csv('compactor_summary_nupack.tsv', sep="\t")
