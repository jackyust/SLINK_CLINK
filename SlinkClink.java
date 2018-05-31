import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

public class SlinkClink {
    public static int algorithmChoice = 1; // 1 for slink, 2 for clink

    public static int numData = 6;
    public static int numDimension = 2;
    public static double cutThreshold = Double.POSITIVE_INFINITY; //Double.POSITIVE_INFINITY;

    public static String directory = "D:\\RSource\\";
    public static String fileName = "slink6objectXY.csv";

    public static double[][] rawData = new double[numData][numDimension];
    public static int[] pi = new int[numData];
    public static double[] lambda = new double[numData];


    public static TreeMap<cluster, ArrayList<Integer>> topLevelCluster = new TreeMap<>(new levelComparator());

    public static void main(String[] args) throws IOException {
        System.out.println("Program Start!");
        long startTime = System.currentTimeMillis();

        readData(directory + fileName);

        if (algorithmChoice == 1)
            slink();
        else if (algorithmChoice == 2)
            clink();

        extractHierarchy();

        long estimatedTime = System.currentTimeMillis() - startTime;
        System.out.println("running time: ");
        System.out.println(estimatedTime);
    }

    public static void slink() throws NumberFormatException, IOException{
        System.out.println("Start Slink algorithm!");

        pi[0] = 1;
        lambda[0] = Double.POSITIVE_INFINITY;
/*		System.out.println("Object 1:");
		System.out.println("pi[1] = " + pi[0]); 
		System.out.println("lambda[1] = " + lambda[0]); 
		System.out.println("================================================================");
*/
        for (int t = 1; t < numData; t++) {
/*			System.out.println("Object:" + (t + 1)); 
			System.out.println("Initialization: pi[" + (t + 1) + "] = " + (t + 1) + ", lambda[" + (t + 1) + "] = " + Double.POSITIVE_INFINITY);
*/			pi[t] = t + 1;
            lambda[t] = Double.POSITIVE_INFINITY;

            //Dynamic temp array for distance between all previous added objects
            double[] M = new double[t];

//			System.out.print("Distance between all previous added objects: ");
            ArrayList<Double> MArrayList = new ArrayList<Double>();
            MArrayList = distanceBetweenPreviousObjects(t);

            for(int MIndex = 0; MIndex < t; MIndex++){
                M[MIndex] = MArrayList.get(MIndex);
            }
			
/*			for (int tempIndex = 0; tempIndex < t; tempIndex++) {
				M[tempIndex] = distanceMatrix[t][tempIndex];				
				System.out.print("M[" + (tempIndex + 1) + "]:" + M[tempIndex]);
				if(tempIndex != t-1)
					System.out.print(", ");
				else
					System.out.println();
			}
*/
            for (int i = 0; i < t; i++) {
//				System.out.println("For previous " + (i + 1) + " objects:");

                if (lambda[i] > M[i]) {
/*					System.out.println("update cases when " + "(lambda[" + (i + 1) + "] = " + lambda[i] + ") > (M[" + (i + 1) + "] = " + M[i] + ")"); 
					System.out.println("M[pi[" + (i + 1) + "]] = Math.min(M[pi[" + (i + 1) + "]], lambda[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.min(M[" + pi[i] + "], lambda[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.min(" + M[pi[i] -1] + ", " + lambda[i] + ") = " + Math.min(M[pi[i] - 1], lambda[i])); 
					System.out.println("M[" + pi[i] + "]: " + M[pi[i] - 1] + " -> " + Math.min(M[pi[i] - 1], lambda[i])); 
*/					M[pi[i] - 1] = Math.min(M[pi[i] - 1], lambda[i]);

//					System.out.println("lambda[" + (i + 1) + "] = M[" + (i + 1) + "] = " + M[i]); 
                    lambda[i] = M[i];

//					System.out.println("pi[" + (i + 1) + "] = " + (t + 1)); 
                    pi[i] = t + 1;
                } else {
/*					System.out.println("update cases when " + "(lambda[" + (i + 1) + "] = " + lambda[i] + ") <= (M[" + (i + 1) + "] = " + M[i] + ")"); 
					System.out.println("M[pi[" + (i + 1) + "]] = Math.min(M[pi[" + (i + 1) + "]], M[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.min(M[" + pi[i] + "], M[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.min(" + M[pi[i] - 1] + ", " + M[i] + ") = " + Math.min(M[pi[i] - 1], M[i])); 
					System.out.println("M[" + pi[i] + "]: " + M[pi[i] - 1] + " -> " + Math.min(M[pi[i] - 1], M[i])); 
*/					M[pi[i] - 1] = Math.min(M[pi[i] - 1], M[i]);
                }
            }

//			System.out.println("Label rearrangement:");
            for (int i = 0; i < t; i++){
                if (lambda[i] >= lambda[pi[i] - 1]){
/*					System.out.println("lambda[" + (i + 1) + "]  = " + lambda[i] + ", lambda[pi[" + (i + 1) + "]] = lambda[" + pi[i] + "] = " + lambda[pi[i] - 1]); //System.out.println("lambda[" + i + "]  = " + lambda[i] + ", lambda[pi[" + i + "]] = lambda[" + pi[i] + "] = " + lambda[pi[i]]);
					System.out.println("pi[" + (i + 1) + "] = " + (t + 1)); 
*/					pi[i] = t + 1;
                }
            }
			
/*			System.out.println("Final status after inserting new object: ");
			for (int i = 0; i <= t; i++){
				System.out.print("pi[" + (i + 1) + "] = " + pi[i]); 
				if(i != t)
					System.out.print(", ");
				else
					System.out.println();			
			}
			
			for (int i = 0; i <= t; i++){								
				System.out.print("lambda[" + (i + 1) + "] = " + lambda[i]); 
				if(i != t)
					System.out.print(", ");
				else
					System.out.println();
			}
			System.out.println("================================================================");
*/		}
		
/*		for (int i = 0; i < numData; i++){
			System.out.println("lambda[" + (i + 1) + "] = " + lambda[i]);
			System.out.println("PI[" + (i + 1) + "] = " + pi[i]);
		}
*/
    }

    public static void clink(){
        System.out.println("Start Clink algorithm!");

        pi[0] = 1;
        lambda[0] = Double.POSITIVE_INFINITY;
/*		System.out.println("Object 1:");
		System.out.println("pi[1] = " + pi[0]); 
		System.out.println("lambda[1] = " + lambda[0]); 
		System.out.println("================================================================");
*/
        for (int t = 1; t < numData; t++) {
/*			System.out.println("Object:" + (t + 1)); 
			System.out.println("Initialization: pi[" + (t + 1) + "] = " + (t + 1) + ", lambda[" + (t + 1) + "] = " + Double.POSITIVE_INFINITY);
*/			pi[t] = t + 1;
            lambda[t] = Double.POSITIVE_INFINITY;

            //Dynamic temp array for distance between all previous added objects
            double[] M = new double[t];

//			System.out.print("Distance between all previous added objects: ");
            ArrayList<Double> MArrayList = new ArrayList<Double>();
            MArrayList = distanceBetweenPreviousObjects(t);

            for(int MIndex = 0; MIndex < t; MIndex++){
                M[MIndex] = MArrayList.get(MIndex);
            }

            for (int i = 0; i < t; i++) {
//				System.out.println("For previous " + (i + 1) + " objects:");

                if(lambda[i] < M[i]){
/*					System.out.println("update cases when " + "(lambda[" + (i + 1) + "] = " + lambda[i] + ") < (M[" + (i + 1) + "] = " + M[i] + ")"); 
					System.out.println("M[pi[" + (i + 1) + "]] = Math.max(M[pi[" + (i + 1) + "]], M[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.max(M[" + pi[i] + "], M[" + (i + 1) + "])"); 
					System.out.println("M[" + pi[i] + "] = Math.max(" + M[pi[i] - 1] + ", " + M[i] + ") = " + Math.max(M[pi[i] - 1], M[i])); 
					System.out.println("M[" + pi[i] + "]: " + M[pi[i] - 1] + " -> " + Math.max(M[pi[i] - 1], M[i]));
					System.out.println("M[" + (i + 1) + "] = " + M[i] + " -> " + "Double.POSITIVE_INFINITY");
*/					M[pi[i] - 1] = Math.max(M[pi[i] - 1], M[i]);
                    M[i] = Double.POSITIVE_INFINITY;
                }
            }

            int a = t - 1;
            for (int i = 0; i < t; i++){
                if (lambda[t - i - 1] >= M[pi[t - i - 1] - 1]){
                    if(M[t - i - 1] < M[a]){
                        a = t - i - 1;
                    }
                }
                else{
                    M[t - i - 1] = Double.POSITIVE_INFINITY;
                }
            }

            int b = pi[a] - 1;
            double c = lambda[a];
            pi[a] = t + 1;
            lambda[a] = M[a];

            if(a < t - 1){
                while(b < t - 1){
                    int d = pi[b] - 1;
                    double e = lambda[b];
                    pi[b] = t + 1;
                    lambda[b] = c;
                    b = d;
                    c = e;
                }

                if(b == t - 1){
                    pi[b] = t + 1;
                    lambda[b] = c;
                }
            }

            for(int i = 0; i < t; i++){
                if(pi[pi[i] - 1] == t + 1 && lambda[i] >= lambda[pi[i] - 1]){
                    pi[i] = t + 1;
                }
            }
			
/*			System.out.println("Final status after inserting new object: ");
			for (int i = 0; i <= t; i++){
				System.out.print("pi[" + (i + 1) + "] = " + pi[i]); 
				if(i != t)
					System.out.print(", ");
				else
					System.out.println();			
			}
			
			for (int i = 0; i <= t; i++){								
				System.out.print("lambda[" + (i + 1) + "] = " + lambda[i]); 
				if(i != t)
					System.out.print(", ");
				else
					System.out.println();
			}
			System.out.println("================================================================");
*/		}
		
/*		for (int i = 0; i < numData; i++){								
			System.out.println("lambda[" + (i + 1) + "] = " + lambda[i]);
			System.out.println("PI[" + (i + 1) + "] = " + pi[i]);
		}
*/
    }

    public static void extractHierarchy() throws IOException{
        System.out.println("extractHierarchical function!");
        int[] orderingIndexByLambda = getOrderingIndex(lambda); //return value: 1-based array(no zero)
        int[] orderedPI = new int[numData];
        double[] orderedLambda = new double[numData];
        orderedPI = intValueOrdering(orderingIndexByLambda, pi); // return value: 1-based(no zero)
        orderedLambda = doubleValueOrdering(orderingIndexByLambda, lambda);

        for(int i = 0; i < numData - 1; i++){
            double height = orderedLambda[i];

            if(height < cutThreshold){
                int level = i + 1;
                cluster c = new cluster(level, height);
                ArrayList<Integer> objectList = new ArrayList<Integer>();

                ArrayList<Integer> blObject = new ArrayList<Integer>();
                ArrayList<Integer> blParent = new ArrayList<Integer>();

                int[] newArray = Arrays.copyOfRange(orderedPI, 0, i);

                int consideredObjectIndex = orderingIndexByLambda[i]; //return value: 1-based
                blObject = checkExistence(consideredObjectIndex, newArray); //return value: may contain 0 value
                blParent = checkExistence(orderedPI[i], newArray);

                int blObjectSize = blObject.size();
                int blParentSize = blParent.size();

                if(blObjectSize > 0 && blParentSize > 0){
                    int objectLastOccurrenceIndex = blObject.get(blObjectSize - 1);
                    objectList.addAll(topLevelCluster.get(new cluster(objectLastOccurrenceIndex + 1, orderedLambda[objectLastOccurrenceIndex])));
                    int piLastOccurrenceIndex = blParent.get(blParentSize - 1);
                    objectList.addAll(topLevelCluster.get(new cluster(piLastOccurrenceIndex + 1, orderedLambda[piLastOccurrenceIndex])));
                }
                else if(blObjectSize > 0 && blParentSize == 0){
                    int objectLastOccurrenceIndex = blObject.get(blObjectSize - 1);
                    objectList.addAll(topLevelCluster.get(new cluster(objectLastOccurrenceIndex + 1, orderedLambda[objectLastOccurrenceIndex])));
                    objectList.add(orderedPI[i]);
                }
                else if(blParentSize > 0 && blObjectSize == 0){
                    int piLastOccurrenceIndex =  blParent.get(blParentSize - 1);
                    objectList.addAll(topLevelCluster.get(new cluster(piLastOccurrenceIndex + 1, orderedLambda[piLastOccurrenceIndex])));
                    objectList.add(consideredObjectIndex);
                }
                else if(blObjectSize == 0 && blParentSize == 0){
                    objectList.add(consideredObjectIndex);
                    objectList.add(orderedPI[i]);
                }
                topLevelCluster.put(c, objectList);
            }
            else
                break;
        }

        Iterator<Map.Entry<cluster, ArrayList<Integer>>> itr = topLevelCluster.entrySet().iterator();

		itr = topLevelCluster.entrySet().iterator();
		  // print all
		while(itr.hasNext()){
			Map.Entry<cluster, ArrayList<Integer>> me = itr.next();
			cluster c = me.getKey();			
			System.out.println("level: " + c.getLevel() + ", height:" + Double.toString(c.getHeight()) + ", list=" + me.getValue());
		}
    }

    public static ArrayList<Double> distanceBetweenPreviousObjects(int consideredObjectIndex){
        ArrayList<Double> list = new ArrayList<Double>();

        for(int i = 0; i < consideredObjectIndex; i++){
            double squaredSum = 0;
            double difference = 0;
            double distance = 0;
            for(int j = 0; j < numDimension; j++){
                difference = rawData[consideredObjectIndex][j] - rawData[i][j];
                squaredSum = squaredSum + Math.pow(difference, 2);
            }
            distance = Math.sqrt(squaredSum);
            list.add(distance);
        }
        return list;
    }

    public static void readData(String filePath) throws NumberFormatException, IOException{
        String line = "";

        // Read distance matrix
        BufferedReader bReader = new BufferedReader(new FileReader(filePath));
        String datavalue[] = null;
        int lineIndex = 0;

        while ((line = bReader.readLine()) != null) {
            datavalue = line.split(",");

            // numData = datavalue.length
            for (int valueIndex = 0; valueIndex < numDimension; valueIndex++)
                rawData[lineIndex][valueIndex] = Double.parseDouble(datavalue[valueIndex]);
            lineIndex++;
        }
        bReader.close();
    }

    public static int[] getOrderingIndex(double[] array) {
        Map<Integer, Double> map = new HashMap<Integer, Double>(array.length);
        for (int i = 0; i < array.length; i++)
            map.put(i, array[i]);

        List<Entry<Integer, Double>> l =
                new ArrayList<Entry<Integer, Double>>(map.entrySet());

        Collections.sort(l, new Comparator<Entry<?, Double>>() {
            @Override
            public int compare(Entry<?, Double> e1, Entry<?, Double> e2) {
                return e1.getValue().compareTo(e2.getValue());
            }
        });

        int[] result = new int[array.length];
        for (int i = 0; i < result.length; i++)
            result[i] = l.get(i).getKey() + 1; //sungho added +1: 20150902

        return result;
    }

    // return value: 0-based
    public static ArrayList<Integer> checkExistence(int targetValue, int arr[]) {
        ArrayList<Integer> indices = new ArrayList<>();
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == targetValue) {
                indices.add(i);
            }
        }
        return indices;
    }

    public static int[] intValueOrdering(int orderIndex[], int arr[]){
        int[] newArray = new int[orderIndex.length];
        for(int i = 0; i < orderIndex.length; i++){
            newArray[i] = arr[orderIndex[i] - 1];
        }
        return newArray;
    }

    public static double[] doubleValueOrdering(int orderIndex[], double arr[]){
        double[] newArray = new double[orderIndex.length];
        for(int i = 0; i < orderIndex.length; i++){
            newArray[i] = arr[orderIndex[i] - 1];
        }
        return newArray;
    }
}

class cluster{
    private int level;
    private double height;

    public cluster(int level, double height) {
        //super();
        this.level = level;
        this.height = height;
    }

    public int getLevel() {
        return level;
    }

    public void setLevel(int level) {
        this.level = level;
    }

    public double getHeight() {
        return height;
    }

    public void setHeight(double height) {
        this.height = height;
    }

    @Override
    public String toString() {
        return "[level=" + level + ", height=" + height + "]";
    }
}

class levelComparator implements Comparator<cluster> {
    @Override
    public int compare(cluster c1, cluster c2) {
        return c2.getLevel() - c1.getLevel(); //descending order
    }
}

class heightComparator implements Comparator<cluster> {
    @Override
    public int compare(cluster c1, cluster c2) {
        if (c1.getHeight() > c2.getHeight()) {
            return 1;
        } else {
            return -1;
        }
    }
}