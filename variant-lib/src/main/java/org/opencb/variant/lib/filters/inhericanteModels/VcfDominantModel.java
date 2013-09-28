package org.opencb.variant.lib.filters.inhericanteModels;


import org.opencb.variant.lib.core.formats.Genotype;
import org.opencb.variant.lib.core.formats.Individual;
import org.opencb.variant.lib.core.formats.Pedigree;
import org.opencb.variant.lib.core.formats.VcfRecord;
import org.opencb.variant.lib.filters.customfilters.VcfFilter;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: aleman
 * Date: 9/27/13
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class VcfDominantModel extends VcfFilter {

    private Pedigree ped;

    public VcfDominantModel(Pedigree ped){
        super();
        this.ped = ped;

    }

    public VcfDominantModel(Pedigree ped, int priority){
        super(priority);
        this.ped = ped;
    }

    @Override
    public boolean apply(VcfRecord vcfRecord) {

        boolean res = false;
        boolean unAff, aff;
        String familyName;
        List<Individual> individuals;

        Iterator<String> it = ped.getFamilies().keySet().iterator();

        while(it.hasNext() && !res){
            familyName = it.next();

            // Check unAffected = 0/0
            individuals = ped.getUnaffected(familyName);
            unAff = true;

            for(Individual ind: individuals){

                Genotype gt = vcfRecord.getSampleGenotype(ind.getId());
                if(gt.isMissing()){
                    ;
                }else if(gt.getAllele1() != 0 || gt.getAllele2() != 0){

                    unAff = false;
                    break;
                }
            }
            if(!unAff){

                System.out.println("Next it");
                continue;
            }

            //   Check Affected = 0/0
            individuals = ped.getAffected(familyName);

            aff = true;

            for(Individual ind: individuals){

                Genotype gt = vcfRecord.getSampleGenotype(ind.getId());
                if(gt.isMissing()){
                    ;
                }
                else if(gt.getAllele1() == 0 && gt.getAllele2() == 0){
                    aff = false;
                    break;
                }

            }

            if(aff){
                res = true;
            }
        }
        return res;
    }
}
