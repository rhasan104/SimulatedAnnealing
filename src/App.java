import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.concurrent.ThreadLocalRandom;
import java.util.Arrays;

class DataSet {
    static int n = 100;
    static int N_t = 1000;
    static int N_m = 100;
    static int N_f = 6 * n;
//    static double c = -1;
//    static double b = 0;
//    static double t = 1.9;
}

class Runner extends Thread {

    private double b;
    private double c;
    private double t;
    public double mean_m = 0;
    public double mean_cp = 0;
    double[] array_m = new double[DataSet.N_m];
    double[] array_Cp = new double[DataSet.N_m];


    Runner(double b, double c, double t){
        this.b = b;
        this.c = c;
        this.t = t;
    }

    public void run(){

        for(int i = 0; i < DataSet.N_m; i++){

            int[] final_config = metropolis_algo();

            array_m[i] = magnetization(final_config);
            array_Cp[i] = pair_correlation(final_config);

        }
        mean_m = getMean(array_m);

        mean_cp = getMean(array_Cp);

    }

    public int[] metropolis_algo(){

        int n = DataSet.n;
        int[] current_config = initial_spin_config();
        int[] new_config = Arrays.copyOf(current_config, n);

        for(int i = 0; i < DataSet.N_f; ++i){
            int s_i = random_spin_position_generator();
            int value_at_s_i = new_config[s_i];

            new_config[s_i] = -value_at_s_i;
            if(update_config(current_config, new_config)) {
                current_config[s_i] = new_config[s_i];
            }
            else {
                new_config[s_i] = value_at_s_i;
            }
        }
        return current_config;
    }

    public int[] initial_spin_config(){

        int n = DataSet.n;
//        double c = DataSet.c;
        int[] config = new int[n];

        if (c >= 0) {
            for(int i = 0; i < n; i++){
                config[i] = 1;
            }
        }
        else {
            int sign = 1;
            for(int i = 0; i < n; i++){
                config[i] = sign;
                sign = -sign;
            }
        }
        return config;
    }

    public int random_spin_position_generator(){
        return ThreadLocalRandom.current().nextInt(DataSet.n);
    }

    public boolean update_config(int[] currentConfig, int[] newConfig){

        double delta_E = 0;
        for(int i = 0; i < DataSet.n; i++){
            int next = i + 1;
            if(next == DataSet.n){
                next = 0;
            }
            delta_E -= (((b * newConfig[i]) - (b * currentConfig[i])) +
                    ((c * newConfig[i] * newConfig[next]) - (c * currentConfig[i] * currentConfig[next])));
        }
        if(delta_E < 0) return true;
        double p = Math.exp(-(delta_E) / t);
        double r = ThreadLocalRandom.current().nextDouble(0, 1);
        if(r < p) return true;

        return false;
    }

    public double magnetization(int[] current_config){
        double magnet = 0;
        for(int i = 0; i < DataSet.n; i++){
            magnet += current_config[i];
        }
        return (1.00 / DataSet.n) * magnet;
    }

    public double pair_correlation(int[] current_config){
        double c_p = 0;
        for(int i = 0; i < DataSet.n; i++){
            int next = i + 1;
            if(next == DataSet.n){
                next = 0;
            }
            c_p += (current_config[i] * current_config[next]);
        }
        return (1.00 / DataSet.n) * c_p;
    }

    public double getMean(double[] value){
        double sum = 0;
        for(int i = 0; i < DataSet.N_m; i++) {
            sum += value[i];
        }
        return (sum / DataSet.N_m);
    }
}

public class App {

    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        BigDecimal bd = new BigDecimal(Double.toString(value));
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }
    public static void findChallenge1() {
        for (double b = 0.1; b <= 2.0; b += round(0.1, 1)) {
            for (double c = -2.0; c <= 0.01; c += round(0.01, 2)) {
                Runner t1 = new Runner(round(b, 1), round(c, 1), round(0.01, 2));
                t1.start();
                try {
                    if (t1.isAlive()) {
                        t1.join();
                    }
                }
                catch (Exception e) {
                    System.out.println("Something went wrong :(");
                }
                if (t1.mean_m == 0 && t1.mean_cp == -1.0) {
                    System.out.println("B = " + b + " & C = " + c + " for:");
                    System.out.println("wantedMeanM = " + t1.mean_m + " & wantedMeanCp = " + t1.mean_cp);
                    return;
                }
            }
        }
        System.out.println("Sorry Couldn't find answer.");
    }

    public static void findChallenge2() {
        for (double b = 0.1; b <= 2.0; b += round(0.1, 1)) {
            for (double c = -2.0; c <= -0.01; c += round(0.01, 2)) {
                Runner t = new Runner(round(b, 1), round(c, 1), round(0.01, 2));
                t.start();
                try {
                    if (t.isAlive()) {
                        t.join();
                    }
                }
                catch (Exception e) {
                    System.out.println("Something went wrong :(");
                }
                if (t.mean_m == 1 && t.mean_cp == 1) {
                    System.out.println("B = " + b + " & C = " + c + " for:");
                    System.out.println("wantedMeanM = " + t.mean_m + " & wantedMeanCp = " + t.mean_cp);
                    return;
                }
            }
        }
        System.out.println("Sorry Couldn't find answer.");
    }

    public static void main(String [] args){

        final long startTime = System.currentTimeMillis();

        int numThreads = DataSet.N_t;
        Runner myT[] = new Runner[numThreads];
        double meanOfM[] = new double[numThreads];
        double meanOfCp[] = new double[numThreads];

        double mu_M = 0;
        double mu_Cp = 0;

        double rel_error = 0;
        double std_div = 0;

//        double cp_theory = ((Math.exp(DataSet.c/DataSet.t) - Math.exp((-DataSet.c)/ DataSet.t))/
//                (Math.exp(DataSet.c/DataSet.t) + Math.exp((-DataSet.c)/DataSet.t)));

        double cp_theory = ((Math.exp(-1/1.9) - Math.exp((1)/ 1.9))/
                (Math.exp(-1/1.9) + Math.exp((1)/1.9)));

        System.out.println("<CP>_theory for b = 0 is " + cp_theory);

        for(int i = 0; i < numThreads; i++){
            myT[i] = new Runner(1.3, -0.05, 47.5);
            myT[i].start();
        }

        for (int i = 0; i < numThreads; i++) {
            try {
                if (myT[i].isAlive()) {
                    myT[i].join();
                }
            }
            catch (Exception e) {
                System.out.println("Something went wrong :(");
            }
        }

        for(int i = 0; i < numThreads; i++){
            meanOfM[i] = myT[i].mean_m;
            meanOfCp[i] = myT[i].mean_cp;
            rel_error += ((myT[i].mean_cp - cp_theory) / cp_theory);

            mu_M += myT[i].mean_m;
            mu_Cp +=myT[i].mean_cp;
        }

        try{
            rel_error = round(rel_error /= numThreads, 4);
            mu_M = round(mu_M /= numThreads, 4);
            mu_Cp = round(mu_Cp /= numThreads, 4);
        } catch (NumberFormatException nfe) {
            System.out.println("NaN dude! try again.");
        }

        for(int i = 0; i < numThreads; i++){
            std_div += Math.pow((((myT[i].mean_cp - cp_theory) / cp_theory) - rel_error), 2);
        }

        try{
            std_div = round(Math.sqrt(std_div /= numThreads), 4);
        } catch (NumberFormatException nfe){
            System.out.println("NaN dude! try again.");
        }

        System.out.println("Global_M = " + mu_M);
        System.out.println("Global_Cp = " + mu_Cp);
        System.out.println("Relative error = " + rel_error);
        System.out.println("Std. Deviation = " + std_div);

//        System.out.println("**************************************************************************");
//        System.out.println("Challenge # 1");
//        findChallenge1();
//
//        System.out.println("**************************************************************************");
//        System.out.println("Challenge # 2");
//        findChallenge2();

        System.out.println("**************************************************************************");

        final long endTime = System.currentTimeMillis();

        System.out.println("Total execution time: " + (endTime - startTime) / 1000.0);

    }
}
