package org.opencb.variant.lib.io;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opencb.variant.lib.core.formats.Family;
import org.opencb.variant.lib.core.formats.Pedigree;
import org.opencb.variant.lib.filters.customfilters.VcfConsequenceTypeFilter;
import org.opencb.variant.lib.filters.customfilters.VcfFilter;
import org.opencb.variant.lib.filters.inhericanteModels.VcfDominantModel;
import org.opencb.variant.lib.io.VariantFilterRunner;
import org.opencb.variant.lib.io.ped.readers.PedDataReader;
import org.opencb.variant.lib.io.ped.readers.PedFileDataReader;
import org.opencb.variant.lib.io.variant.annotators.VcfAnnotator;
import org.opencb.variant.lib.io.variant.annotators.VcfConsequenceTypeAnnotator;
import org.opencb.variant.lib.io.variant.annotators.VcfGeneNameAnnotator;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: javi
 * Date: 26/09/13
 * Time: 18:44
 * To change this template use File | Settings | File Templates.
 */
public class VariantFilterRunnerTest {


    @Before
    public void setUp() throws Exception {



    }

    @After
    public void tearDown() throws Exception {

    }

//    @Test
//    public void filterCT(){
//        String fileIn = "/home/javi/tmp/small_annot.vcf";
//        String fileOut = "/home/javi/tmp/small_out.vcf";
//
//        VariantFilterRunner vr = new VariantFilterRunner(fileIn, fileOut);
//        List<VcfFilter> list = new ArrayList<>();
//
//
//        list.add(new VcfConsequenceTypeFilter("DNAseI_hypersensitive_site"));
//
//
//        vr.filters(list);
//
//        vr.run();
//
//
//    }

    @Test
    public void filterDominant(){
        String fileIn = "/opt/data/mgp/dominant/dominant.vcf";
        String fileOut = "/opt/data/mgp/dominant/dominant_out.vcf";
        String pedFile = "/opt/data/mgp/dominant/dominant.ped";

        PedDataReader pedReader = new PedFileDataReader(pedFile);

            pedReader.open();
            Pedigree ped = pedReader.read();
            pedReader.close();

        VariantFilterRunner vr = new VariantFilterRunner(fileIn, fileOut);

        List<VcfFilter> list = new ArrayList<>();


        list.add(new VcfDominantModel(ped));


        vr.filters(list);

        vr.run();

    }
}
