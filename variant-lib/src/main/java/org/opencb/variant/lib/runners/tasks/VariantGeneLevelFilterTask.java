package org.opencb.variant.lib.runners.tasks;

import org.opencb.commons.bioformats.pedigree.Condition;
import org.opencb.commons.bioformats.pedigree.Individual;
import org.opencb.commons.bioformats.pedigree.Pedigree;
import org.opencb.commons.bioformats.pedigree.io.readers.PedigreePedReader;
import org.opencb.commons.bioformats.variant.Variant;
import org.opencb.commons.bioformats.variant.VariantStudy;
import org.opencb.commons.bioformats.variant.filters.VariantCompoundHeterozygosityFilter;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecordCoordinateComparator;
import org.opencb.commons.bioformats.variant.filters.VariantGeneLevelFilter;
import org.opencb.commons.bioformats.variant.vcf4.io.writers.VariantWriter;
import org.opencb.commons.run.Task;
import org.opencb.variant.lib.runners.VariantRunner;

import java.io.IOException;
import java.util.*;

/**
 * Created by parce on 2/27/14.
 */
public class VariantGeneLevelFilterTask extends Task<Variant>{


    private Pedigree pedigree;

    private List<VariantGeneLevelFilter> filters;

    private VariantWriter writer;

    private static final int UNAFFECTED = 1;

    private static final int AFFECTED = 2;

    // TODO: este campo esta ya en el runner, estudiar si sobra aqui
    private VariantStudy study;

    public VariantGeneLevelFilterTask(VariantStudy study, List<VariantGeneLevelFilter> geneFilters, VariantWriter writer) {
        this(study, geneFilters, writer, null);
    }

    public VariantGeneLevelFilterTask(VariantStudy study, List<VariantGeneLevelFilter> geneFilters, VariantWriter writer, VariantRunner prev)
    {
        //super(study, reader, pedReader, writer, prev);
        super();
        this.study = study;
        this.pedigree = this.study.getPedigree();
        this.filters = geneFilters;
        this.writer = writer;
    }

    public boolean apply(List<Variant> batch) throws IOException {
        List<Variant> tempBatch;
        for (VariantGeneLevelFilter filter : this.filters) {
            tempBatch = filter.apply(batch);
            batch = tempBatch;
        }
        // TODO: hacer que el metodo apply de los filtros filtre directamente sobre el batch de entrada
        return true;
    }

//    // NUEVA VERSION DEL METODO:
//    // Usamos el reader que va devolviendo batchs de variantes por gen, asi que el batch que recibe este metodo ya son
//    // las variantes de un unico gen, y podemos pasar directamente al filtro
//    public boolean apply(List<Variant> batch) throws IOException {
//        //List<Variant> filteredBatch = new LinkedList<Variant>();
//        // TODO: cambiar nombre
//        Set<Variant> preOutputSet = new TreeSet<>(new VcfRecordCoordinateComparator());
//        // TODO: Hay VcfRecordCoordinateComparator de la clase Variant?
//        List<Variant> variantsToRemove = null;
//        List<String> finishedGenes;
//        List<String> geneNames;
//        for (Variant variant : batch) {
//
//                // filter the variants from finished genes, adding those who pass the filters to the output batch
//                preOutputSet.addAll(this.applyGeneLevelFilterToFinishedGenes(finishedGenes, geneMap));
//
//                // add to the output batch (in coordinate order) the variants that are not in "not finished" genes
//                boolean found = false;
//                for (Variant record : preOutputSet) {
//                    for (List<Variant> geneVariants: geneMap.values()) {
//                        if (geneVariants.contains(record)) {
//                            found = true;
//                            break;
//                        }
//                    }
//                    if (!found) {
//                        filteredBatch.add(record);
//                        if (variantsToRemove == null) {
//                            variantsToRemove = new LinkedList<Variant>();
//                        }
//                        variantsToRemove.add(record);
//                    } else {
//                        break;
//                    }
//                }
//                // remove the processed variants from the set
//                if (variantsToRemove != null) {
//                    preOutputSet.removeAll(variantsToRemove);
//                    variantsToRemove = null;
//                }
//
//        }
//        // write the filtered batch
//        batch.clear();
//        if (writer != null) {
//            ((VariantWriter) writer).write(filteredBatch);
//        }
//
//        // TODO: ver si es necesario devolver una variable (true si va bien, false si hay error)
//        return true;
//
//    }

//    public boolean apply(List<Variant> batch) throws IOException {
//        List<Variant> filteredBatch = new LinkedList<Variant>();
//        // TODO: cambiar nombre
//        Set<Variant> preOutputSet = new TreeSet<Variant>(new VcfRecordCoordinateComparator());
//        // TODO: Hay VcfRecordCoordinateComparator de la clase Variant?
//        List<Variant> variantsToRemove = null;
//        List<String> finishedGenes;
//        List<String> geneNames;
//        for (Variant variant : batch) {
//            // check if there are "finished" genes
//            geneNames = this.variantGeneNames(variant);
//            finishedGenes = this.finishedGenes(geneNames, geneMap.keySet());
//            if (finishedGenes != null) {
//                // filter the variants from finished genes, adding those who pass the filters to the output batch
//                preOutputSet.addAll(this.applyGeneLevelFilterToFinishedGenes(finishedGenes, geneMap));
//                // remove the finished genes from the gene map
//                for (String finishedGene : finishedGenes) {
//                    geneMap.remove(finishedGene);
//                }
//                // add to the output batch (in coordinate order) the variants that are not in "not finished" genes
//                boolean found = false;
//                for (Variant record : preOutputSet) {
//                    for (List<Variant> geneVariants: geneMap.values()) {
//                        if (geneVariants.contains(record)) {
//                            found = true;
//                            break;
//                        }
//                    }
//                    if (!found) {
//                        filteredBatch.add(record);
//                        if (variantsToRemove == null) {
//                            variantsToRemove = new LinkedList<Variant>();
//                        }
//                        variantsToRemove.add(record);
//                    } else {
//                        break;
//                    }
//                }
//                // remove the processed variants from the set
//                if (variantsToRemove != null) {
//                    preOutputSet.removeAll(variantsToRemove);
//                    variantsToRemove = null;
//                }
//            }
//
//            // if the variant has one of more gene names, add them to the gene map
//            if (geneNames != null) {
//                this.addVariantAndGenes(variant, geneNames);
//            }
//        }
//        // write the filtered batch
//        batch.clear();
//        if (writer != null) {
//            ((VariantWriter) writer).write(filteredBatch);
//        }
//
//        // TODO: el prototipo dice que debe devolver un boolean. Co
//        return filteredBatch;
//    }

//
//    /**
//     *
//     * @param finishedGenes
//     * @param geneMap
//     * @return
//     */
//    private Set<Variant> applyGeneLevelFilterToFinishedGenes(List<Variant> variants) {
//        Set<Variant> passedVariants = new HashSet<Variant>();
//        //Map<String, Individual> individuals = this.pedigree.getIndividuals();
//        for (String gene : finishedGenes) {
//            passedVariants.addAll(this.filter.apply(geneMap.get(gene)));
//        }
//        return passedVariants;
//    }
//
//    /**
//     * Add the variant and its genes to the geneMap
//     * @param variant - variant to be added
//     * @param geneNames - List with one or more Gene names
//     */
//    private void addVariantAndGenes(Variant variant, List<String> geneNames) {
//        List<Variant> variantsFromGene;
//        for (String geneName : geneNames) {
//            variantsFromGene = geneMap.get(geneName);
//            if (variantsFromGene == null) {
//                // the gene is not in the map, create a new variant list
//                variantsFromGene = new LinkedList<Variant> ();
//            }
//            variantsFromGene.add(variant);
//            geneMap.put(geneName, variantsFromGene);
//        }
//    }
//
//    // TODO: este metodo deberia estar en el PED
//    private Set<String> getIndividualsWithCondition(Condition condition) {
//        Set<String> filteredIndividualNames = new HashSet<String>();
//        Map<String, Individual> individuals = this.pedigree.getIndividuals();
//
//        for (String individualName : individuals.keySet()) {
//            if (individuals.get(individualName).getCondition() == condition) {
//                filteredIndividualNames.add(individualName);
//            }
//        }
//
//        return filteredIndividualNames;
//    }
}
