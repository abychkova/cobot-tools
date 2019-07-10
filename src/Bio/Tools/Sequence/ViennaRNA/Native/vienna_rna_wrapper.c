#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>

float
vrna_fold_temperature(
  const char *string,
  char *structure,
  double temperature) {

    float                 mfe;
    vrna_fold_compound_t  *vc;
    vrna_md_t             md;

    vrna_md_set_default(&md);
    md.temperature = temperature;
    vc  = vrna_fold_compound(string, &md, 0);
    mfe = vrna_mfe(vc, structure);

    vrna_fold_compound_free(vc);

    return mfe;
}

float
vrna_cofold_temperature(
  const char *seq,
  char *structure,
  double temperature) {

    float                 mfe;
    vrna_fold_compound_t  *vc;
    vrna_md_t             md;

    vrna_md_set_default(&md);
    md.temperature = temperature;
    md.min_loop_size = 0;  /* set min loop length to 0 */

    /* get compound structure */
    vc = vrna_fold_compound(seq, &md, 0);

    mfe = vrna_mfe_dimer(vc, structure);

    vrna_fold_compound_free(vc);

    return mfe;
}
