
import java.io.*;
import java.util.*;

public class CHUI_PSO {

    private Map<Integer, Integer> itemTWU = new HashMap<>(); //map of item TWUs
    private List<List<Pair>> database = new ArrayList<>(); //the pruned database
    private Particle gBest; //global fittest particle
    private Particle[] pBest; //list of fittest particle offspring
    private Particle[] population; //the particle population
    private int maxTransactionLength = 0; //longest transaction length in database
    private ArrayList<Item> items = new ArrayList<>(); //the 1-HTWUI
    private List<Particle> chuis = new ArrayList<>(); //solution-set
    private HashSet<BitSet> explored = new HashSet<>(); //set of explored particles
    private HashMap<Integer, Integer> itemNamesRev = new HashMap<>(); //reference to old original item names
    private int shortestTransactionID; //The shortest transaction of the current particle
    private int std; //the mean deviation between avg- and max utilities
    private boolean avgEstimate; //boolean that decides the type of estimate
    private int lowEst = 0; //underestimates
    private int highEst = 0; //overestimates



    //file paths
    final String dataset = "retail";
    final String input = "D:\\Documents\\Skole\\Master\\Work\\" + dataset + ".txt"; //input file path
    final String output = "D:\\Documents\\Skole\\Master\\Work\\out.txt"; //output file path

    //Algorithm parameters
    final int pop_size = 20; // the size of the population
    final int iterations = 10000; // the number of iterations before termination
    final int minUtil = 5000; // minimum utility threshold
    final boolean closed = true; //true = find CHUIS, false = find HUIS
    final boolean prune = true; //true = ETP, false = traditional TWU-Model

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
        BitSet TIDS; //TidSet
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
     * Call this method to run the algorithm. File paths and algorithm parameters must be set in top of this class
     *
     * @throws IOException
     */
    public void run() throws IOException {
        maxMemory = 0;
        startTimestamp = System.currentTimeMillis();

        readData(); //reads input file and prunes DB
        checkMemory();


        //utilities used after each population update
        List<Double> probChui = new ArrayList<>(); //roulette probabilities for current discovered CHUIs
        int nPatterns = 0; // number of CHUIs discovered last iteration
        int pos = 0; //position of gBest in percentChui
        int lastImproved = 0; //the number of iterations since a CHUI was discovered
        boolean roulette = true; // roulette wheel selection

        //calculate average utility of each item and find the deviation between avgUtil & maxUtil
        std = 0; // the standard deviation
        for (Item item : items) {
            item.avgUtil = 1 + (item.totalUtil / item.TIDS.cardinality());
            std += item.maxUtil - item.avgUtil;
        }

        explored.add(new BitSet(items.size())); //avoids edge-case for empty particle
        if (!items.isEmpty()) {
            std = std / items.size(); //find mean deviation
            //only use avgEstimates if the standard deviation is small compared to the minUtil
            avgEstimate = (double) std / minUtil < 0.0001;
            generatePop(); //initialize the population
            for (int i = 0; i < iterations; i++) {
                //update and evaluate each particle in population
                update();

                //gBest update strategy
                if (chuis.size() > 1) {
                    //SELECT GBEST WITH RWS
                    if (roulette) {
                        if (nPatterns != chuis.size()) { //new CHUIs discovered, roulette probs must be updated
                            probChui = rouletteProbChui();
                            lastImproved = 0;
                        } else {
                            lastImproved++;
                            if (lastImproved > 50) {
                                roulette = false; //starting to converge, disable RWS
                            }
                        }
                        pos = rouletteSelect(probChui); //select gBest with RWS
                    }
                    //SELECT GBEST AS NEXT CHUI IN 'chuis'
                    else {
                        if (nPatterns == chuis.size()) { //only change gBest if no CHUIs discovered last iteration
                            if (pos < chuis.size() - 1) {
                                pos++;
                            } else {
                                pos = 0; //end of list reached, start at front again
                            }
                        }
                    }
                    gBest = new Particle(chuis.get(pos).X, chuis.get(pos).fitness); //update gBest
                }


                if (i % 100 == 0 && highEst > 0 && i > 0) { //check each 100th iteration
                    //Tighten std if mostly overestimates are made (only relevant when avgEstimates is active)
                    std = ((double) lowEst / highEst < 0.01) ? 1 : std;
                }
                nPatterns = chuis.size();
            }
        }
        endTimestamp = System.currentTimeMillis();
        checkMemory();
        writeOut();
    }


    /**
     * Generates pop_size particles with roulette wheel selection based on item TWUs
     */
    private void generatePop() {
        List<Double> rouletteProbabilities = rouletteProbabilities();
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
            BitSet tidSet = pev_check(p); //transactions the particle occur
            p.fitness = calcFitness(p, tidSet, -1);
            population[i] = p;
            pBest[i] = new Particle(p.X, p.fitness); //initialize pBest
            if (!explored.contains(p.X)) {
                //check if HUI/CHUI
                if (p.fitness >= minUtil) {
                    if (closed) {
                        //check if particle is closed
                        if (isClosed(p, tidSet)) {
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

    /**
     * The pev-check verifies that the particle exists in the database and modifies it if not.
     * Furthermore, it calculates the avg/max fitness estimate and returns the TidSet of the particle.
     *
     * @param p The particle
     * @return orgBitSet: The transactions the itemset of the particle occur (TidSet)
     */
    private BitSet pev_check(Particle p) {
        int item1 = p.X.nextSetBit(0);
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
                // no common transactions, remove the current item from the particle
                orgBitSet = (BitSet) copyBitSet.clone();
                p.X.clear(i);
            }
        }
        return orgBitSet;
    }

    /**
     * Verifies the closure of a particle
     *
     * @param p                   The particle
     * @param tidSet              the TIDSET of the particle
     * @return True if Closed, false otherwise
     */
    private boolean isClosed(Particle p, BitSet tidSet) {
        int support = tidSet.cardinality();
        List<Pair> sTrans = database.get(shortestTransactionID);
        //Loop all items in the shortest transaction the particle occur
        for (Pair pair : sTrans) {
            //The item does not appear in the particle
            if (!p.X.get(pair.item)) {
                BitSet currentTids = (BitSet) tidSet.clone();
                BitSet newItem = (BitSet) items.get(pair.item - 1).TIDS.clone();
                currentTids.and(newItem);
                //The support of the particle is the same with the new item appended -> The particle is not closed
                if (currentTids.cardinality() == support) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Calculates the fitness of a particle
     *
     * @param p      The particle
     * @param tidSet TidSet of the particle
     * @param idx    The position of the particle in the population (to reference pBest), set to -1 if first population
     * @return The fitness of the particle
     */
    private int calcFitness(Particle p, BitSet tidSet, int idx) {
        int fitness = 0;
        int min = maxTransactionLength; //used to find the shortest transaction

        //The particle only contains 1 item, return the fitness calculated during pre-processing
        if (p.X.cardinality() == 1) {
            return items.get(p.X.nextSetBit(0) - 1).totalUtil;
        }

        //estimate the fitness
        int support = tidSet.cardinality();
        int est = p.estFitness * support;
        int buffer = avgEstimate ? (std * support) : 0;
        if (idx != -1) {
            if (est + buffer < minUtil && est < pBest[idx].fitness) {
                // Skip fitness calculation
                return 0;
            }
        }

        //calculate exact fitness
        for (int i = tidSet.nextSetBit(0); i != -1; i = tidSet.nextSetBit(i + 1)) {
            int q = 0; //current index in transaction
            int item = p.X.nextSetBit(0);
            if (database.get(i).size() < min) {
                shortestTransactionID = i; //update shortest transaction id of particle
            }
            while (item != -1) {
                if (database.get(i).get(q).item == item) {
                    fitness += database.get(i).get(q).utility;
                    item = p.X.nextSetBit(item + 1);
                }
                q++;
            }
        }

        //Update overestimates and underestimates
        if (est + buffer < fitness) {
            lowEst++;
        } else {
            highEst++;
        }
        return fitness;
    }


    /**
     * Updates population and evaluates the population
     */
    private void update() {
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
                BitSet tidSet = pev_check(p);

                //check if explored again because pev_check can change the particle
                if (!explored.contains(p.X)) {
                    p.fitness = calcFitness(p, tidSet, i);
                    //update pBest and gBest
                    if (p.fitness > pBest[i].fitness) {
                        pBest[i] = new Particle(p.X, p.fitness);
                        if (p.fitness > gBest.fitness) {
                            gBest = new Particle(p.X, p.fitness);
                        }
                    }
                    if (p.fitness >= minUtil) {
                        if (closed) {
                            if (isClosed(p, tidSet)) {
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
     * Flips a random number of bits in current particle, only bits that are opposite to pBest/gBest are considered
     *
     * @param diffList bit differences between particle and pBest/gBest
     * @param pos      position of current particle in population
     */
    private void changeParticle(List<Integer> diffList, int pos) {
        //number of items so change
        int num = (int) (diffList.size() * Math.random() + 1);
        if (diffList.size() > 0) {
            for (int i = 0; i < num; i++) {
                //position to change
                int change = (int) (diffList.size() * Math.random());
                //flip the bit of the selected item and remove it from difflist
                population[pos].X.flip(diffList.get(change));
            }
        }
    }

    /**
     * Computes the bit difference between current particle and pBest/gBest
     *
     * @param best pBest/gBest
     * @param p    the current particle
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
     * @return List of probability ranges
     */
    private List<Double> rouletteProbabilities() {
        List<Double> probRange = new ArrayList<>();
        double twuSum = 0;
        double tempSum = 0;
        //sum the twu values for all 1-HTWUIs
        for (int i = 0; i < items.size(); i++) {
            int item = items.get(i).item;
            twuSum += itemTWU.get(item);
        }
        //Set probabilities based on TWU-proportion
        for (int i = 0; i < items.size(); i++) {
            int item = items.get(i).item;
            tempSum += itemTWU.get(item);
            double percent = tempSum / twuSum;
            probRange.add(percent);
        }
        return probRange;
    }

    /**
     * Select an element based on the given roulette wheel probabilities
     *
     * @param probRange list of roulette probabilities
     * @return position of selected element in list
     */
    private int rouletteSelect(List<Double> probRange) {
        double rand = Math.random();
        if (rand <= probRange.get(0)) {
            return 0;
        }
        int pos = 0;
        for (int i = 1; i < probRange.size(); i++) {
            if (rand > probRange.get(i - 1) && rand <= probRange.get(i)) {
                pos = i;
                break;
            }
        }
        return pos;
    }

    /**
     * creates roulette probabilities from current CHUIs
     * @return list of probabilities for each CHUI
     */
    private List<Double> rouletteProbChui() {
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
     * Reads DB from input file and initializes the pruning strategies + transaction optimizations
     */
    private void readData() {
        Map<Integer, Integer> itemTWU1 = new HashMap<>(); //holds current TWU-value for each item
        List<Integer> transUtils = new ArrayList<>(); //holds TU-value for each transaction
        List<List<Pair>> db = new ArrayList<>();
        List<List<Pair>> tempDb = new ArrayList<>();
        String currentLine;
        try (BufferedReader data = new BufferedReader(new InputStreamReader(
                new FileInputStream(input)))) {
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
            ETP(db, transUtils); //Use additional pruning with ETP
        } else {
            optimizeTransactions(db, itemTWU1);
        }

    }

    /**
     * Recursively calculates item-TWUs, removes 1-LTWUI and updates TUs, until no items are removed.
     *
     * @param db         The database to prune
     * @param transUtils The current transaction utilities of the database
     */
    private void ETP(List<List<Pair>> db, List<Integer> transUtils) {
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
                } else { // item is 1-LTWUI
                    pruned = true;
                    int TU = transUtils.get(i);
                    TU -= db.get(i).get(j).utility;
                    transUtils.set(i, TU); //update transaction utility
                }
            }
            revisedDB.add(revisedTransaction); //store the revised transaction
        }
        if (pruned) { //item was removed, repeat pruning
            ETP(revisedDB, transUtils);
        } else { //pruning is finished, optimize DB
            optimizeTransactions(revisedDB, itemTWU1);
        }
    }

    /**
     * Sets item-names in the range 1 - #1-HTWUI, removes empty transactions,
     * and initializes values required for the fitness- calculation and estimate approach
     * The revised db is stored in 'database'
     *
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
     * Writes results to file output file
     * @throws IOException
     */
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
        BufferedWriter w = new BufferedWriter(new FileWriter(output));
        w.write(sb.toString());
        w.newLine();
        w.close();
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

    /**
     * check current memory usage
     */
    private void checkMemory() {
        double currentMemory = (Runtime.getRuntime().totalMemory() - Runtime
                .getRuntime().freeMemory()) / 1024d / 1024d;
        maxMemory = Math.max(maxMemory, currentMemory);
    }
}
