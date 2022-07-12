
import java.io.*;
import java.util.*;

public class CHUI_PSO {

    private Map<Integer, Integer> itemTWU = new HashMap<>();
    private List<List<Pair>> database = new ArrayList<>();
    private Particle gBest;
    private Particle[] pBest;
    private Particle[] population;
    private int maxTransactionLength = 0;
    private ArrayList<Item> items = new ArrayList<>();
    private List<Particle> chuis = new ArrayList<>();
    private HashSet<BitSet> explored = new HashSet<>();
    private HashMap<Integer, Integer> itemNamesRev = new HashMap<>();
    private int shortestTransactionID;
    private int avgDiff;
    private boolean avgEstimate;
    private int lowEst = 0;
    private int highEst = 0;
    ArrayList<Integer> it = new ArrayList<>();
    ArrayList<Integer> pat = new ArrayList<>();



    //file paths
    final String dataset = "accidents";
    final String dataPath = "D:\\Documents\\Skole\\Master\\Work\\"+dataset+".txt"; //input file path
    final String resultPath = "D:\\Documents\\Skole\\Master\\Work\\out.txt"; //output file path
    final String convPath = "D:\\Documents\\Skole\\Master\\Experiments\\"+dataset+"\\";

    //Algorithm parameters
    final int pop_size = 20; // the size of the population
    final int iterations = 10000; // the number of iterations
    final int minUtil = 29000000; // minimum utility threshold
    final boolean closed = true; //true = find CHUIS, false = find HUIS
    final boolean prune = true; //true = extended pruning, false = traditional pruning

    //stats
    double maxMemory = 0; // the maximum memory usage
    long startTimestamp = 0; // the time the algorithm started
    long endTimestamp = 0; // the time the algorithm terminated
    long totalUtil = 0; // the total utility in the DB


    // this class represent an item and its utility in a transaction
    private static class Pair implements Comparable<Pair> {
        int item;
        int utility;

        public Pair(int item, int utility) {
            this.item = item;
            this.utility = utility;
        }

        @Override
        public int compareTo(Pair o) {
            return (this.item < o.item) ? -1 : 1;
        }
    }

    private static class Item {
        int item; //item-name
        BitSet TIDS; //TID set
        int totalUtil = 0; //utility of the item
        int avgUtil; // average utility of item
        int maxUtil = 0; // maximum utility of item

        public Item(int item) {
            TIDS = new BitSet();
            this.item = item;
        }

        @Override
        public String toString() {
            return String.valueOf(item);
        }
    }

    // this class represent the particles
    private static class Particle {
        BitSet X; // items contained in particle (encoding vector)
        int fitness; // fitness/utility of particle
        int estFitness; // estimated fitness of particle

        public Particle(int size) {
            this.X = new BitSet(size);
        }

        public Particle(BitSet bitset, int fitness) {
            this.X = (BitSet) bitset.clone();
            this.fitness = fitness;
        }

        @Override
        public String toString() {
            return X.toString();
        }

    }


    /**
     * @throws IOException
     */
    public void run() throws IOException {
        maxMemory = 0;
        startTimestamp = System.currentTimeMillis();

        readData(); //reads input file and prunes DB

        System.out.println("db size: " + database.size());
        System.out.println("twu size: " + items.size());
        System.out.println("trans: " + maxTransactionLength);
        System.out.println(totalUtil);

        checkMemory();

        //variables used after each population update
        List<Double> percentChui = new ArrayList<>(); //roulette probabilities for current CHUIs
        int nPatterns = 0; // number of patterns discovered last iteration
        int pos = 0; //position of gBest in percentChui
        int lastImproved = 0; //the number of iterations since a CHUI was discovered
        boolean roulette = true; // roulette wheel selection

        //calculate average utility of each item and find the average difference between them all
        avgDiff = 0;
        for (Item item : items) {
            item.avgUtil = 1 + (item.totalUtil / item.TIDS.cardinality());
            avgDiff += item.maxUtil - item.avgUtil;
        }

        if (!items.isEmpty()) {
            avgDiff = avgDiff / items.size();
            System.out.println("avg_diff: " + avgDiff);
            //if the proportion of the average difference to the minUtil is too large,
            // then use maximum utilities for fitness estimates
            avgEstimate = (double) avgDiff / minUtil < 0.0001;
            if(avgEstimate) {
                System.out.println("avg_est");
            }
            else {
                System.out.println("max_est");
            }
            generatePop(); //initialize the population
            for (int i = 0; i < iterations; i++) {
                update(); //update each particle in population
                //gBest update strategy
                if (chuis.size() > 1) {
                    if (roulette) {
                        if (nPatterns != chuis.size()) { //new CHUIs discovered
                            percentChui = roulettePercentChui(); //recalculate roulette probabilities
                            lastImproved = 0;
                        } else {
                            lastImproved++;
                            if (lastImproved > 50) {
                                roulette = false; //convergence has slowed, disable roulette wheel selection
                            }
                        }
                        pos = rouletteSelect(percentChui); //select gBest
                    } else { //select gBest with different approach
                        if (nPatterns == chuis.size()) { //only change gBest if no CHUIs discovered last iteration
                            if (pos < chuis.size() - 1) {
                                pos++;
                            } else {
                                pos = 0;
                            }
                        }
                    }
                    gBest = new Particle(chuis.get(pos).X, chuis.get(pos).fitness); //update gBest
                }
                if(nPatterns != chuis.size()) {
                    it.add(i);
                    pat.add(chuis.size());
                    System.out.println("iteration: " + i + " CHUIs: " + chuis.size());
                }

                if (i % 100 == 0 && highEst > 0) { //check each 100th iteration
                    //Tighten avgDiff if mostly overestimates are made (only relevant when avgEstimates is active)
                    avgDiff = ((double) lowEst / highEst < 0.01) ? 1 : avgDiff;
                }
                nPatterns = chuis.size();
            }
        }
        endTimestamp = System.currentTimeMillis();
        checkMemory();
        writeOut();
        //writeRes();
    }


    private void generatePop() { //TODO: can be more effective when creating already explored particles.
        List<Double> rouletteProbabilities = roulettePercent();
        population = new Particle[pop_size];
        pBest = new Particle[pop_size];
        for (int i = 0; i < pop_size; i++) {
            //k is the number of items to include in the particle
            int k = (int) (Math.random() * maxTransactionLength) + 1;
            //j is the current number of items that has been included
            int j = 0;
            Particle p = new Particle(items.size());
            while (j < k) {
                //select the item
                int pos = rouletteSelect(rouletteProbabilities);
                //if item is not previously selected, select it and increment j.
                if (!p.X.get(items.get(pos).item)) {
                    p.X.set(items.get(pos).item);
                    j++;
                }
            }
            BitSet transactionBitset = pev_check(p); //transactions the particle occur
            p.fitness = calcFitness(p, transactionBitset, -1);
            population[i] = p;
            pBest[i] = new Particle(p.X, p.fitness); //initialize pBest
            if (!explored.contains(p.X)) {
                //check if HUI/CHUI
                if (p.fitness >= minUtil) {
                    if (closed) {
                        if (isClosed(p, shortestTransactionID, transactionBitset)) {
                            chuis.add(new Particle(p.X, p.fitness));
                        }
                    } else {
                        chuis.add(new Particle(p.X, p.fitness));
                    }
                }
            }
            if (i == 0) {
                gBest = new Particle(p.X, p.fitness); //initialize gBest
            } else {
                if (p.fitness > gBest.fitness) {
                    gBest = new Particle(p.X, p.fitness); //update gBest
                }
            }
            BitSet clone = (BitSet) p.X.clone();
            explored.add(clone); //set particle as explored
        }
    }

    private BitSet pev_check(Particle p) {
        int item1 = p.X.nextSetBit(0);
        if (item1 == -1) {
            return null; //TODO: create new particle?? and why does this seem to always happen once
        }
        BitSet orgBitSet = (BitSet) items.get(item1 - 1).TIDS.clone();
        BitSet copyBitSet = (BitSet) orgBitSet.clone();
        p.estFitness = avgEstimate ? items.get(item1 - 1).avgUtil : items.get(item1 - 1).maxUtil;
        for (int i = p.X.nextSetBit(item1 + 1); i != -1; i = p.X.nextSetBit(i + 1)) {
            orgBitSet.and(items.get(i - 1).TIDS);
            //the two items have common transactions
            if (orgBitSet.cardinality() > 0) {
                copyBitSet = (BitSet) orgBitSet.clone();
                p.estFitness += avgEstimate ? (items.get(i - 1).avgUtil) : (items.get(i - 1).maxUtil);
            } else {
                // no common transactions
                orgBitSet = (BitSet) copyBitSet.clone();
                p.X.clear(i);
            }
        }
        return orgBitSet;
    }

    private boolean isClosed(Particle p, int shortestTransaction, BitSet Tids) {
        int support = Tids.cardinality();
        List<Pair> lTrans = database.get(shortestTransaction);
        for (Pair pair : lTrans) {
            if (!p.X.get(pair.item)) {
                BitSet currentTids = (BitSet) Tids.clone();
                BitSet newItem = (BitSet) items.get(pair.item - 1).TIDS.clone();
                currentTids.and(newItem);
                if (currentTids.cardinality() == support) {
                    return false;
                }
            }
        }
        return true;
    }

    private int calcFitness(Particle p, BitSet transactionBitset, int idx) {
        int fitness = 0;
        int min = maxTransactionLength;
        if (transactionBitset == null) {
            return 0;
        }

        if (p.X.cardinality() == 1) {
            return items.get(p.X.nextSetBit(0) - 1).totalUtil;
        }
        //TODO: rewrite
        //estimate the fitness and determine if is worth calculating exact fitness
        int support = transactionBitset.cardinality();
        int est = p.estFitness * support;
        //delta is used to account for wobble in the estimate
        double delta = minUtil - (avgDiff * support);
        delta = (delta < 0) ? 0 : delta;
        //only use delta if avgEstimate is used. Using max utility values ensures est >= fitness
        delta = avgEstimate ? (delta / minUtil) : 1;
        if (idx != -1) {
            //don't calculate fitness if estimate is unpromising
            if (est < minUtil * delta && est < pBest[idx].fitness) {
                return 0;
            }
        }

        //calculate exact fitness
        for (int i = transactionBitset.nextSetBit(0); i != -1; i = transactionBitset.nextSetBit(i + 1)) {
            int q = 0; //current index in transaction
            int item = p.X.nextSetBit(0);
            if (database.get(i).size() < min) {
                shortestTransactionID = i; //TODO: Better to store a set with each item that does not occur in P
            }
            while (item != -1) {
                if (database.get(i).get(q).item == item) {
                    fitness += database.get(i).get(q).utility;
                    item = p.X.nextSetBit(item + 1);
                }
                q++;
            }
        }

        if (est < fitness) {
            lowEst++;
        } else {
            highEst++;
        }
        return fitness;
    }


    /**
     * Updates population and checks for new CHUIs
     */
    private void update() { //TODO: Memory opt.
        List<Integer> diffList;
        for (int i = 0; i < pop_size; i++) {
            Particle p = population[i];
            //different bits between pBest and current particle
            diffList = bitDiff(pBest[i], p);
            //change a random amount of these bits
            changeParticle(diffList, i);
            //repeat for gBest
            diffList = bitDiff(gBest, p);
            changeParticle(diffList, i);


            if (explored.contains(p.X)) {
                //the particle is already explored, change one random bit
                int rand = (int) (items.size() * Math.random());
                int change = items.get(rand).item;
                p.X.flip(change);
            }

            //avoid PEV-check and fit. calc. if particle is already explored
            if (!explored.contains(p.X)) {
                //bitset before pev
                BitSet copy1 = (BitSet) p.X.clone();
                BitSet transactionBitset = pev_check(p);

                //check if explored again because pev_check can change the particle
                if (!explored.contains(p.X)) {
                    p.fitness = calcFitness(p, transactionBitset, i);
                    //update pBest and gBest
                    if (p.fitness > pBest[i].fitness) {
                        pBest[i] = new Particle(p.X, p.fitness);
                        if (p.fitness > gBest.fitness) {
                            gBest = new Particle(p.X, p.fitness);
                        }
                    }
                    if (p.fitness >= minUtil) {
                        if (closed) {
                            if (isClosed(p, shortestTransactionID, transactionBitset)) {
                                //Particle is CHUI
                                chuis.add(new Particle(p.X, p.fitness));
                            }
                        } else {
                            // particle is HUI
                            chuis.add(new Particle(p.X, p.fitness));
                        }
                    }
                    //bitset after pev
                    BitSet copy2 = (BitSet) p.X.clone();
                    explored.add(copy2); //set current particle as explored
                }
                explored.add(copy1); //set particle before PEV-check as explored
            }
        }
    }

    /**
     *  Flips a random number of bits in current particle, only bits that a different to pBest/gBest are considered
     * @param diffList bit differences between particle and pBest/gBest
     * @param pos position of current particle in population
     */
    private void changeParticle(List<Integer> diffList, int pos) {
        //number of items so change
        int num = (int) (diffList.size() * Math.random() + 1);
        if (diffList.size() > 0) {
            for (int i = 0; i < num; i++) {
                //position to change
                int change = (int) (diffList.size() * Math.random());
                //flip the bit of the selected item
                population[pos].X.flip(diffList.get(change));
            }
        }
    }

    /**
     * Computes the bit difference between current particle and pBest/gBest
     * @param best pBest/gBest
     * @param p the current particle
     * @return List containing the different bit positions
     */
    private List<Integer> bitDiff(Particle best, Particle p) {
        List<Integer> diffList = new ArrayList<>();
        BitSet temp = (BitSet) best.X.clone();
        temp.xor(p.X);
        for (int i = temp.nextSetBit(0); i != -1; i = temp.nextSetBit(i + 1)) {
            diffList.add(i);
        }
        return diffList;
    }



    /**
     * Creates a list of probabilities for roulette wheel selection based on item TWU-values
     *
     * @return
     */
    private List<Double> roulettePercent() {
        List<Double> percents = new ArrayList<>();
        double twuSum = 0;
        double tempSum = 0;
        //sum the twu values for all 1-HTWUIs
        for (int i = 0; i < items.size(); i++) {
            int item = items.get(i).item;
            twuSum += itemTWU.get(item);
        }
        //Calculate the proportional probabilities
        for (int i = 0; i < items.size(); i++) {
            int item = items.get(i).item;
            tempSum += itemTWU.get(item);
            double percent = tempSum / twuSum;
            percents.add(percent);
        }
        return percents;
    }

    /**
     * Select the
     *
     * @param percents list of roulette percentages
     * @return
     */
    private int rouletteSelect(List<Double> percents) {
        double rand = Math.random();
        if (rand <= percents.get(0)) {
            return 0;
        }
        int pos = 0;
        for (int i = 1; i < percents.size(); i++) {
            if (rand > percents.get(i - 1) && rand <= percents.get(i)) {
                pos = i;
                break;
            }
        }
        return pos;
    }

    private List<Double> roulettePercentChui() {
        double sum = 0;
        double tempSum = 0;
        List<Double> percentsChui = new ArrayList<>();
        for (Particle hui : chuis) {
            sum += hui.fitness;
        }
        for (Particle hui : chuis) {
            tempSum += hui.fitness;
            double percent = tempSum / sum;
            percentsChui.add(percent);
        }
        return percentsChui;
    }


    /**
     * Reads transactions from input file and removes items that are not 1-HTWUIs
     * if prune=true: extended pruning will be initialized
     * if prune=false: DB is revised, optimized and returned in 'database' list
     */
    private void readData() {
        Map<Integer, Integer> itemTWU1 = new HashMap<>(); //holds current TWU-value for each item
        List<Integer> transUtils = new ArrayList<>(); //holds TU-value for each transaction
        List<List<Pair>> db = new ArrayList<>();
        List<List<Pair>> tempDb = new ArrayList<>();
        String currentLine;
        try (BufferedReader data = new BufferedReader(new InputStreamReader(
                new FileInputStream(dataPath)))) {
            //1st DB-Scan: calculate TWU value for each item
            while ((currentLine = data.readLine()) != null) {
                String[] split = currentLine.split(":");
                String[] items = split[0].split(" ");
                String[] utilities = split[2].split(" ");
                int transactionUtility = Integer.parseInt(split[1]);
                totalUtil += transactionUtility;
                transUtils.add(transactionUtility);
                List<Pair> transaction = new ArrayList<>();
                for (int i = 0; i < items.length; i++) {
                    int item = Integer.parseInt(items[i]);
                    int util = Integer.parseInt(utilities[i]);
                    Pair pair = new Pair(item, util);
                    transaction.add(pair);
                    Integer twu = itemTWU1.get(item);
                    twu = (twu == null) ? transactionUtility : twu + transactionUtility;
                    itemTWU1.put(item, twu);
                }
                tempDb.add(transaction);
            }
        } catch (Exception e) {
            // catches exception if error while reading the input file
            e.printStackTrace();
        }
        //2nd DB-scan: remove items with TWU < minUtil
        for (int i = 0; i < tempDb.size(); i++) {
            List<Pair> revisedTransaction = new ArrayList<>();
            for (int j = 0; j < tempDb.get(i).size(); j++) {
                Pair pair = tempDb.get(i).get(j);
                if (itemTWU1.get(pair.item) >= minUtil) {
                    revisedTransaction.add(pair);
                } else {
                    int TU = transUtils.get(i);
                    TU -= pair.utility;
                    transUtils.set(i, TU); //update transaction utility since item is removed
                }
            }
            db.add(revisedTransaction); //store revised transaction
        }
        if (prune) {
            prune(db, transUtils); //extended pruning
        } else {
            optimizeTransactions(db, itemTWU1);
        }

    }

    /**
     * Sets item-names in the range 1 - #candidate-items, removes empty transactions,
     * and initializes values required for the fitness- calculation and estimate approach
     * Can significantly reduce memory usage and execution time with certain datasets
     * @param db The database to optimize
     */
    private void optimizeTransactions(List<List<Pair>> db, Map<Integer, Integer> itemTWU1) {
        HashMap<Integer, Integer> itemNames = new HashMap<>();
        int c = 0; //new item name
        int transID = 0; //current TID
        for (int i = 0; i < db.size(); i++) {
            if (db.get(i).isEmpty()) {
                continue;
            }
            for (int j = 0; j < db.get(i).size(); j++) {
                int item = db.get(i).get(j).item;
                int utility = db.get(i).get(j).utility;
                if (!itemNames.containsKey(item)) {
                    //item has not been given new name yet
                    c++; //increment name
                    itemNames.put(item, c); //set name for this item
                    itemNamesRev.put(c, item); //save the old name so it can be retrieved later
                    Item itemClass = new Item(c); //this class stores different info for the item
                    items.add(itemClass);
                }
                int twu = itemTWU1.get(item); //get the twu of this item
                item = itemNames.get(item); //get the name of the item
                db.get(i).get(j).item = item; //change the name of the item
                Item it = items.get(item - 1);
                it.TIDS.set(transID); //update the transaction bit for the item
                itemTWU.put(item, twu); //store twu value
                it.totalUtil += utility; //update total utility of this item
                it.maxUtil = (it.maxUtil == 0) ? utility : Math.max(it.maxUtil, utility); //update max utility
            }
            Collections.sort(db.get(i)); //sort transaction according to item name
            maxTransactionLength = Math.max(maxTransactionLength, db.get(i).size()); //update max trans. length
            database.add(db.get(i)); //store the transaction
            transID++;
        }
    }

    /**
     * Recursively calculates item-TWUs, prune unpromising items, and update TUs, until no items are removed.
     * @param db         The database to prune
     * @param transUtils The current transaction utilities of the database
     */
    private void prune(List<List<Pair>> db, List<Integer> transUtils) {//TODO: make iterative

        List<List<Pair>> revisedDB = new ArrayList<>();
        Map<Integer, Integer> itemTWU1 = new HashMap<>();
        boolean pruned = false;
        //calculate TWU of each item
        for (int i = 0; i < db.size(); i++) {
            int transactionUtility = transUtils.get(i);
            for (int j = 0; j < db.get(i).size(); j++) {
                int item = db.get(i).get(j).item;
                Integer twu = itemTWU1.get(item);
                twu = (twu == null) ? transactionUtility : twu + transactionUtility;
                itemTWU1.put(item, twu);
            }
        }
        //check if any item has TWU < minUtil
        for (int i = 0; i < db.size(); i++) {
            List<Pair> revisedTransaction = new ArrayList<>();
            for (int j = 0; j < db.get(i).size(); j++) {
                int item = db.get(i).get(j).item;
                int twu = itemTWU1.get(item);
                if (twu >= minUtil) {
                    revisedTransaction.add(db.get(i).get(j)); //add item to revised transaction
                } else { // item is not 1-HTWUI
                    pruned = true;
                    int TU = transUtils.get(i);
                    TU -= db.get(i).get(j).utility;
                    transUtils.set(i, TU); //update transaction utility
                }
            }
            revisedDB.add(revisedTransaction); //store the revised transaction
        }
        if (pruned) { //item was removed, prune again
            prune(revisedDB, transUtils);
        }
        else { //pruning is finished, optimize DB
            optimizeTransactions(revisedDB, itemTWU1);
        }
    }


    private void writeOut() throws IOException {
        StringBuilder sb = new StringBuilder();
        for (Particle p : chuis) {
            for (int i = p.X.nextSetBit(0); i != -1; i = p.X.nextSetBit(i + 1)) {
                sb.append(itemNamesRev.get(i));
                sb.append(" ");
            }
            sb.append("#UTIL: ");
            sb.append(p.fitness);
            sb.append(System.lineSeparator());
        }
        BufferedWriter w = new BufferedWriter(new FileWriter(resultPath));
        w.write(sb.toString());
        w.newLine();
        w.close();
    }

    private void writeRes() throws IOException {
        StringBuilder sb = new StringBuilder();
        String p;
        if(prune) {
            p = "_PRUNE";
        }
        else {
            p = "_NOPRUNE";
        }
        for(int i = 0; i < it.size(); i++) {
            sb.append(it.get(i));
            sb.append(",");
            sb.append(pat.get(i));
            sb.append(System.lineSeparator());
        }
        BufferedWriter w = new BufferedWriter(new FileWriter(convPath+"convergence"+p+".csv"));
        w.write(sb.toString());
        w.newLine();
        w.close();

        sb = new StringBuilder();
        sb.append(minUtil + "," + (minUtil * 1.0 / totalUtil)  + "," + chuis.size() +"," +(endTimestamp - startTimestamp)+ "," + (int) maxMemory);
        BufferedWriter s = new BufferedWriter(new FileWriter(convPath+"log_"+p+".csv", true));
        s.write(sb.toString());
        s.newLine();
        s.close();
    }

    /**
     * Print statistics about the latest execution to System.out.
     */
    public void printStats() {
        System.out
                .println("============= STATS =============");
        System.out.println(" Total time ~ " + (endTimestamp - startTimestamp)
                + " ms");
        System.out.println(" Memory ~ " + maxMemory + " MB");
        System.out.println(" Closed High-utility itemsets count : " + chuis.size());
        System.out
                .println("===================================================");
    }

    private void checkMemory() {
        double currentMemory = (Runtime.getRuntime().totalMemory() - Runtime
                .getRuntime().freeMemory()) / 1024d / 1024d;
        maxMemory = Math.max(maxMemory, currentMemory);
    }
}
