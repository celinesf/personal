
import java.io.*;
import com.mongodb.*;
import java.util.*; 
import java.lang.reflect.Method;  
import java.net.*;  
import java.util.jar.*;  
import java.util.zip.GZIPInputStream;


///*** Value class to have pointer on hapmap data ***/
//class HashMap<String,HashMap<String,Object>> {
//	HashMap<String,HashMap<String,Object>> data ;
//}
///*** Value class to have pointer on haplotype data ***/
//class HashMap<String,List<HashMap<String,Object>>> {
//	HashMap<String,List<HashMap<String,Object>>> data ;
//}

///*** Value class to have pointer on haplotypeInfo data ***/
//class HashMap<String,HashMap<String,Object>> {
//	HashMap<String,HashMap<String,Object>> data ;
//}

/***
 * 
 * 
 * @author celine
 *java -Xmx991m -jar beagle.jar
  missing=?
  markers=/Users/celine/genophen/genophyzer/UnphasedData/markers_chr10
  out=/Users/celine/genophen/genophyzer/PhasedData/phased
  unphased=/Users/celine/genophen/genophyzer/UnphasedData/unphased_chr10
  trios=/Users/celine/genophen/genophyzer/UnphasedData/trios0_chr10
  trios=/Users/celine/genophen/genophyzer/UnphasedData/trios1_chr10
  trios=/Users/celine/genophen/genophyzer/UnphasedData/trios2_chr10
  unphased=/Users/celine/genophen/genophyzer/UnphasedData/unrelated0_chr10
 */

public class Haplophyzer {
	public static void main(String[] args) throws Exception {

        /* initiation for pre-processing */
        String homeDir = System.getProperty("user.dir");//TODO change the path to input/output data
        String inputFile = "output.txt";
        String fileUnphased="unphased";
		String dirUnphased="UnphasedData";
		String filePhased="phased";
		String dirPhased="PhasedData";
		String phasedStart="phased.unphased_chr";
		String phasedEnd="phased.gz";
		String outputFile = "haplotype" ;
		String trios="trios";
		String marker = "markers_chr";
		String unrelated="unrelated";
		HashMap<String,List<HashMap<String,Object>>> halotypeList = new HashMap<String,List<HashMap<String,Object>>>();
		HashMap<String,HashMap<Double,Object>> snpPosition = new HashMap<String,HashMap<Double,Object>>();// chr/position hasmap
		HashMap<String,List<Double>> chromosomePosition = new HashMap<String,List<Double>>();// chr/ sorted list of positions
        /* end of init */
	    System.out.println("I am running haplophyzer on file dir:"+ homeDir + "/"+ inputFile);
	    
	    /* pre processing */
	    Haplophyzer.prepareBeableInput(homeDir,inputFile,dirUnphased,fileUnphased, halotypeList,snpPosition,chromosomePosition) ;
	    System.out.println(halotypeList);
	    
	    /* calling beagle with object as input */
	    Haplophyzer.runBeagle(homeDir,marker,trios,unrelated,dirUnphased,dirPhased, filePhased);
	    
	    /* process beagle output */
	    Haplophyzer.processBeagleOutput(homeDir, inputFile, outputFile, dirPhased,phasedStart, phasedEnd,halotypeList);

	    System.out.println("---> Haplophyzer is done for: " +inputFile);
    }/*** END MAIN ***/
	
	
	/*** prepareBeableInput ***/
	private static void prepareBeableInput(String dirData,String inputFile,String dirUnphased,String fileUnphased,HashMap<String,List<HashMap<String,Object>>> halotypeList,HashMap<String, HashMap<Double,Object>> snpPosition,HashMap<String,List<Double>> chromosomePosition){
		System.out.println("---> prepareBeableInput for Phasing file: " +inputFile);

		/*** inialization ***/
		List<HashMap<String,Object>> trioData = new ArrayList<HashMap<String,Object>>() ;
		List<HashMap<String,Object>> unrelatedData = new ArrayList<HashMap<String,Object>>();

		/*** get HapMap data ***/
		getHapmapData(trioData,unrelatedData);
		
		/*** get list of haplotypes and snps data ***/
		getHaplotypeInfo(halotypeList, snpPosition,chromosomePosition);
		
		/*** get the member genotype for these SNPs ***/
		getMemberGenotype( dirData, inputFile, snpPosition);
		System.out.println(snpPosition);
		
		/*** create input files for beagle ***/
		File dir = new File(dirData +"/"+ dirUnphased);  
		dir.mkdir();
		for (String chromosome :  chromosomePosition.keySet() ){
			generateBeagleInputFiles(chromosome, chromosomePosition.get(chromosome), snpPosition.get(chromosome),trioData,unrelatedData, dirData, dirUnphased, fileUnphased);
		}	
	}/*** END prepareBeableInput ***/

	/*** getMemberGenotype ***/
	private static void getMemberGenotype(String dirData,String inputFile, HashMap<String,HashMap<Double,Object>> snpPosition){
		System.out.println("---> getMemberGenotype ");
		/*** get list of SNPs + positions ***/
		HashMap<String,HashMap<String,Double>> snpMap = new HashMap<String,HashMap<String,Double>>();
		for (String chromosome: snpPosition.keySet()){
			snpMap.put(chromosome, new HashMap<String,Double>());
			for (Double position: snpPosition.get(chromosome).keySet()){
				@SuppressWarnings("unchecked")
				HashMap<String,String> dbsnpInfo = (HashMap<String,String> ) snpPosition.get(chromosome).get(position);
				String  dbsnp = dbsnpInfo.get("dbsnp");
				snpMap.get(chromosome).put(dbsnp, position);
			}
		}
		System.out.println("here\n" + snpMap); 
		
		try { 
    	    /* open file */
    		FileInputStream input = new FileInputStream( dirData+"/"+ inputFile);	   
    	    InputStreamReader xover = new InputStreamReader(input);
			BufferedReader fileBuffer = new BufferedReader(xover);
			
			/* read file / get genotype */
			String line = "";
			while (( line = fileBuffer.readLine()) != null) {
				String[] word = line.split("\t");
				String chromosome = word[0];
				if(snpPosition.containsKey(chromosome) && snpMap.containsKey(chromosome)){
					String dbsnp = (word[2]);
					if(snpMap.get(chromosome).containsKey(dbsnp)){
						String genotype = word[3];
						Double position = snpMap.get(chromosome).get(dbsnp);
						@SuppressWarnings("unchecked")
						HashMap<String,String> snpInfo = (HashMap<String,String>) snpPosition.get(chromosome).get(position);
						if (snpInfo.get("dbsnp").equals(dbsnp)){
							snpInfo.put("gt", genotype);
							snpPosition.get(chromosome).put(position,snpInfo);
						}
						else{
							System.out.println("PROBLEM" + chromosome + " " + dbsnp + " but at this position I find " +snpInfo.get("dbsnp") );
						}						
					}/* if dbsnp in list */
				}/* if chromosome in lists */
			}/* loop read file */
			fileBuffer.close();
        } 
        catch (Exception e) {
        	System.out.println("NO FILE " + dirData+"/"+ inputFile);
			throw new RuntimeException();
		}
	}/*** END prepareBeableInput ***/
	
	private static void generateBeagleInputFiles(String chromosome,List<Double> chromosomePosition,  HashMap<Double,Object> snpPosition, List<HashMap<String,Object>>  trioData,List<HashMap<String,Object>>  unrelatedData , String dirData,String dirUnphased,String fileUnphased){
		System.out.println("---> generateBeagleInputFiles for chromosome " + chromosome);
		
		/*** init ***/
		File markerFile = new File(dirData + "/" +  dirUnphased + "/markers_chr" + chromosome);
		File beagleFile = new File( dirData + "/" +  dirUnphased + "/" + fileUnphased + "_chr" + chromosome);
		
		/*** get/create list of files with hapmap data ***/
		HashMap<String,List<HashMap<String,Object>>>  trioFile = new HashMap<String,List<HashMap<String,Object>>> ();
		createHapmapInputFile(chromosome,trioFile,  trioData,  dirData,dirUnphased, "trios", 0,0);
		HashMap<String,List<HashMap<String,Object>>>  unrelatedFile = new HashMap<String,List<HashMap<String,Object>>> ();
		createHapmapInputFile(chromosome,unrelatedFile,  unrelatedData,  dirData,dirUnphased, "unrelated", 0,0);
		try {
			/* marker */
			FileWriter markerInput = new FileWriter(markerFile.getAbsoluteFile());
			BufferedWriter markerWriter = new BufferedWriter(markerInput);
			/* unphase member data */
			
			FileWriter beagleInput = new FileWriter(beagleFile.getAbsoluteFile());
			BufferedWriter beagleWriter = new BufferedWriter(beagleInput);
			beagleWriter.write("I\tid\tmember\tmember\n")   ;
			
			for (Double position : chromosomePosition ){
				@SuppressWarnings("unchecked")
				HashMap<String,String> snpInfo= (HashMap<String,String>) snpPosition.get(position);
				String genotype = snpInfo.get("gt");
				if (snpInfo.get("gt").equals("NA")){
					genotype = "??"; 
				}
				String dbsnp = snpInfo.get("dbsnp");
				String majorAllele = snpInfo.get("major_allele");
				String minorAllele = snpInfo.get("minor_allele");			
				markerWriter.write(dbsnp+"\t"+chromosome+"\t"+ (Long)position.longValue()+"\t"+majorAllele +"\t"+minorAllele+"\n");
				beagleWriter.write("M\t" + dbsnp +" \t" + genotype.charAt(0) +" \t"+ genotype.charAt(1) +" \n");
				fillHapmapInputFile(snpInfo, trioFile);
				fillHapmapInputFile(snpInfo, unrelatedFile);
	        }
			beagleWriter.close();
			markerWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}/*** END generateBeagleInputFiles ***/

	
	/*** fillHamapInputFile ***/
	private static void fillHapmapInputFile(HashMap<String,String> snpInfo, HashMap<String,List<HashMap<String,Object>>>  fileList){
		System.out.println("---> fillHamapInputFile ");
		for (String fileName :fileList.keySet()){
			File outputFile = new File(fileName);
			try {
				/* marker */
				FileWriter output = new FileWriter(outputFile.getAbsoluteFile(),true);
				BufferedWriter writer = new BufferedWriter(output);
				/** write header **/
				writer.write("M\t" +snpInfo.get("dbsnp") +"\t");		
				for (HashMap<String,Object> snpData : fileList.get(fileName)){
					@SuppressWarnings("unchecked")
					HashMap<String,String> snpList = (HashMap<String,String> )snpData.get("genotype");
					String dbsnp = snpInfo.get("dbsnp");
					String genotype = snpList.get(dbsnp);
					if (genotype != null && !genotype.equals("NN")){
						writer.write(genotype.charAt(0) +" \t"+ genotype.charAt(1) +" \t");
						if (dbsnp.equals("rs7903146")){
						System.out.println(dbsnp+" "+ genotype.charAt(0)+" "+ genotype.charAt(1)+" "+ snpData.get("indv1"));
						}
	                }
                    else{
                    	writer.write("?\t?\t");
	                }
				}
				writer.write("\n");
				writer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} 
	}/*** END fillHamapInputFile ***/
	
	
	/*** createHapmapInputFile ***/
	private static void createHapmapInputFile(String chromosome,HashMap<String,List<HashMap<String,Object>>> fileList, List<HashMap<String,Object>> hapmapData,String dirData,String dirUnphased, String nameFile, Integer nunFile,Integer numIndv){
		System.out.println("---> createHapmapInputFile for "+ nameFile +" "+nunFile+" "+ numIndv);
		
		String outputName = dirData+"/"+ dirUnphased +"/"+ nameFile +nunFile +"_chr"+ chromosome;
		fileList.put(outputName, new ArrayList<HashMap<String,Object>>());

		File outputFile = new File(outputName);
		try {
			/* marker */
			FileWriter output = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter writer = new BufferedWriter(output);
			/** write header **/
			writer.write("I\tid\t")   ;			
			Integer i = 0;
			for (i = numIndv; (i < numIndv+60 && i < hapmapData.size()) ; i++){
				String indv = hapmapData.get(i).get("indv1").toString();
				fileList.get(outputName).add( hapmapData.get(i));
				writer.write(indv +"\t"+indv+"\t" );
			}
			if(i == numIndv+60){
				writer.write("\n");
				writer.close();
				nunFile ++;
				createHapmapInputFile(chromosome,fileList,  hapmapData,  dirData,dirUnphased, nameFile, nunFile,i);
			}
			else{
				writer.write("\n");
				writer.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}/*** END createHapmapInputFile ***/
	    		
	    		
	/*** getHapmapData ***/
	private static void getHapmapData(List<HashMap<String,Object>> trioData,List<HashMap<String,Object>> unrelatedData){
		System.out.println("---> getHapmapData  ");
		HashMap<String,HashMap<String,Object>>  hapmapData = getHapmapIndvData();
		
		getHapmapTrio(trioData,hapmapData) ;

		getHapmapUnrelated(unrelatedData,hapmapData);
	}/*** END getHapmapData ***/

	
	/*** getHapmapIndvData ***/
	private static  HashMap<String,HashMap<String,Object>>  getHapmapIndvData(){
		System.out.println("---> getHapmapIndvData  ");
		HashMap<String,HashMap<String,Object>> hapmap = new HashMap<String,HashMap<String,Object>>();	
		try {
			/*** Connect to Mongodb TODO ***/ 
			Mongo mongo = new Mongo("50.112.142.17", 27017);
			DB db = mongo.getDB("genetics");
			db.authenticate("adm1n", "adm1npw0rd".toCharArray());
			
			/*** Recover hapmap individual names and information ***/
			DBCollection collectionIndv = db.getCollection("HapMapFamily");
			DBCursor HapMapFamily = collectionIndv.find();
			
			while (HapMapFamily.hasNext()) {
				DBObject dbObject = HapMapFamily.next();
				String indv1 = dbObject.get("indv1").toString();
				@SuppressWarnings (value="unchecked")
				HashMap<String,Object> doc = (HashMap<String,Object>) dbObject.toMap();
				hapmap.put(indv1,doc);
			}
			/*** Recover hapmap genotypeing data ***/
			DBCollection collectionData = db.getCollection("FamilyAllele");
			DBCursor FamilyAllele = collectionData.find();
			
			while (FamilyAllele.hasNext()) {
				DBObject dbObject = FamilyAllele.next();
				String indv1 = dbObject.get("Indv").toString();
				if (!hapmap.get(indv1).containsKey("genotype")) {
					HashMap<String,String> genotype = new HashMap<String,String>();
					genotype.put(dbObject.get("dbSNP").toString(), dbObject.get("A1").toString()+dbObject.get("A2").toString());	
					hapmap.get(indv1).put("genotype",genotype);
				}
				else{
					@SuppressWarnings (value="unchecked")
					HashMap<String,String> genotype = (HashMap<String,String>) hapmap.get(indv1).get("genotype");
					genotype.put(dbObject.get("dbSNP").toString(), dbObject.get("A1").toString()+dbObject.get("A2").toString());	
					hapmap.get(indv1).put("genotype",genotype);
				}
			}
			 
		} catch (UnknownHostException e) {
			e.printStackTrace();
		} catch (MongoException e) {
			e.printStackTrace();
		}
		return hapmap;
	}/*** END getHapmapIndvData ***/
	
	
	/*** getHapmapTrio ***/
	private static void getHapmapTrio(List<HashMap<String,Object>> trioData,HashMap<String,HashMap<String,Object>> hapmap){
		System.out.println("---> getHapmapTrio  ");

		for (String hapmapIndv : hapmap.keySet()) {
			if (hapmap.get(hapmapIndv).containsKey("genotype")) {
				String mother = hapmap.get(hapmapIndv).get("indv3").toString() ;
				String father = hapmap.get(hapmapIndv).get("indv2").toString() ;				
				if (!father.equals("0")  &&  !mother.equals("0")){
					getMumAndDad(hapmapIndv,mother,father,trioData,hapmap);
	            }
				/*** I have dad info ***/
            	else if (!father.equals("0") ) {
            		getParentOnly(hapmapIndv,father,hapmap,1);
            	}
            	else if (!mother.equals("0")){ 
            		getParentOnly(hapmapIndv,mother,hapmap,2);
            	}
            	else{
            		/*** do nothing -> no mum or dad***/
            		}
			}
			else{
				/*** individual don't have any genotype data ***/
				System.out.println("WARNING: I have no genotype data for individual " + hapmapIndv);
			}
		}/* loop on hapmap individual */
	}/*** END getHapmapTrio ***/
	
	
	/*** getParentOnly ***/
	private static void	getParentOnly(String hapmapIndv,String parent, HashMap<String,HashMap<String,Object>> hapmap,Integer gender){
//		System.out.println("---> getParentOnly  ");
		String relation = "dad";
		if (gender == 2){
			relation = "mul";
		}
   		if (( Integer ) hapmap.get(parent).get("Gender") == gender 
    			&& hapmap.get(parent).get("Family").toString().equals(hapmap.get(hapmapIndv).get("Family").toString())){   
			/** only dad data **/
        	if (hapmap.get(parent).get("genotype") != null) {
        		hapmap.get(hapmapIndv).put("in", 0);   // only use the father unrelated 
            }/** should not be here? **/
        	else {
        		System.out.println("WARNING: no data for " + relation + " of " + hapmap.get(parent));
        	}
		}
		else{
			System.out.println("ISSUE with Family " + hapmap.get(hapmapIndv).get("Family") + " of indv " +hapmapIndv+ "'s " + relation + " is "+ parent);
		}
	}/*** END getParentOnly ***/
	
	
	/*** getMumAndDad ***/
	private static void	getMumAndDad(String hapmapIndv,String mother,String father,List<HashMap<String,Object>> trioData,HashMap<String,HashMap<String,Object>> hapmap){
//		System.out.println("---> getMumAndDad  ");
		if (( Integer ) hapmap.get(father).get("Gender") == 1 
    			&& ( Integer )hapmap.get(mother).get("Gender") == 2 
    			&& hapmap.get(mother).get("Family").toString().equals(hapmap.get(hapmapIndv).get("Family").toString())
    			&& hapmap.get(father).get("Family").toString().equals(hapmap.get(hapmapIndv).get("Family").toString())){   
	    	/** add a full new trio **/
			if (hapmap.get(mother).get("genotype") != null && hapmap.get(father).get("genotype") != null) {   
				hapmap.get(hapmapIndv).put("in", 0);
	    		hapmap.get(mother).put("in", 0);
	    		hapmap.get(father).put("in", 0);
	    		trioData.add(hapmap.get(mother));
	    		trioData.add(hapmap.get(father));
	    		trioData.add(hapmap.get(hapmapIndv));
	    	}
			/** only mum data **/
	    	else if (hapmap.get(mother).get("genotype") != null) {
	    		hapmap.get(hapmapIndv).put("in", 0);
	    	}
			/** only dad data**/
	    	else if (hapmap.get(father).get("genotype") != null) {
	    		hapmap.get(hapmapIndv).put("in", 0);   
	        }
			/** should not be here? **/
	    	else {
	    		System.out.println("PROBLEM key " + hapmapIndv + " family " + hapmap.get(hapmapIndv).get("Family") + " dad " + hapmap.get(father) + " mum " +  hapmap.get(mother));
	    	}
		}
		/** should not be here? **/
    	else{
    	  System.out.println("ISSUE with Family " + hapmap.get(hapmapIndv).get("Family") + " of indv " +hapmapIndv);
    	}
	}/*** END getMumAndDad ***/
		
	
	/*** getHapmapUnrelated ***/
	private static void  getHapmapUnrelated(List<HashMap<String,Object>> unrelatedData, HashMap<String,HashMap<String,Object>> hapmap){
		System.out.println("---> getHapmapUnrelated  ");
		for (String key : hapmap.keySet()) {
			if (hapmap.get(key).containsKey("genotype") && !hapmap.get(key).containsKey("in")) {
				hapmap.get(key).put("in",0);
				unrelatedData.add(hapmap.get(key));
			}
		}
	}
	/*** END getHapmapUnrelated ***/
	
	
	/*** runBeagle ***/
	private static void runBeagle(String homeDir,String marker,String trios,String unrelated,String dirUnphased,String dirPhased,String filePhased) throws Exception {
		System.out.println("---> runBeagle ");
		
		/** loop on chromosome numbers **/
		(new File(homeDir+"/"+ dirPhased)).mkdirs();
		File[] listUnphasedFiles = new File(homeDir + "/"+dirUnphased).listFiles();
		for (File unphasedFile : listUnphasedFiles) {
		    if (unphasedFile.isFile() && unphasedFile.getName().startsWith("unphased_chr")) {
		    	/* prepare cmd line */
		    	String chrNum =  unphasedFile.getName().split("_chr")[1];
		    	System.out.println(unphasedFile.getName() + " "+ chrNum);
		        String[] cmdLine= new String[0];
		        cmdLine = addElement(cmdLine, "-Xmx991m" + "" + "beagle.jar"); // TODO ADD BEAGLE DIRECTORY
		        cmdLine = addElement(cmdLine,"missing=?");
		        cmdLine = addElement(cmdLine,"markers=" + homeDir +"/"+ dirUnphased+"/"+ marker + chrNum);
		        cmdLine = addElement(cmdLine,"out=" + homeDir +"/"+ dirPhased+"/"+ filePhased) ; 
		        cmdLine = addElement(cmdLine,"unphased=" + homeDir +"/"+ dirUnphased+"/"+ unphasedFile.getName());
		        cmdLine = addHapmapFiles(homeDir,chrNum,cmdLine, dirUnphased,trios,"trios");
		        cmdLine = addHapmapFiles(homeDir,chrNum,cmdLine, dirUnphased,unrelated,"unphased");
		    	
				/* run beagle TODO change path to beagle from conf file*/
			    File jarfile = new File("beagle.jar");  
			    JarFile jar = new JarFile(jarfile);  
			    Manifest manifest = jar.getManifest();  
			    Attributes attrs = manifest.getMainAttributes();  
			    String mainClassName = attrs.getValue("Main-Class");  
			    URL url = new URL("file", null, jarfile.getAbsolutePath());  
			    ClassLoader cl = new URLClassLoader(new URL[] {url});  
		    	Class<?> mainClass = cl.loadClass(mainClassName);  
			    Method mainMethod = mainClass.getMethod("main", new Class[] {String[].class});  
			    String[] arguments = new String[cmdLine.length - 1];  
			    System.arraycopy(cmdLine, 1, arguments, 0, cmdLine.length - 1);  
			    mainMethod.invoke(mainClass, new Object[] {arguments});
		    }
		}
		deleteDirectory(new File(homeDir +"/"+ dirUnphased));
	}/*** END runBeagle ***/
	
	/*** addElement 
	 * Add string at the end of an array
	 * ***/
	private static String[] addElement(String[] org, String added) {
		System.out.println("---> addElement added:" + added);
		String[] result = Arrays.copyOf(org, org.length +1);
	    result[org.length] = added;
	    return result;
	}/*** END addElement ***/

	/*** addHapmapFiles 
	 * Add arguments related to hapmap for beagle run
	 * ***/
	private static String[] addHapmapFiles(String homeDir,String chrNum,String[] cmdLine,String dirData,String typeFile,String typeBeagle){
		System.out.println("---> addHapmapFiles ");
		File[] listFiles = new File(homeDir + "/"+ dirData).listFiles();
		for (File file : listFiles) {
		    if (file.isFile() && file.getName().startsWith(typeFile) && file.getName().endsWith("_chr" +chrNum)) {
		    	cmdLine = addElement(cmdLine,typeBeagle +"="+ homeDir +"/"+ dirData+"/"+ file.getName());
		    }
		}
		return cmdLine;
	}/*** END addElement ***/
	
	
	/*** processBeagleOutput ***/
	private static void processBeagleOutput(String homeDir,String inputFile,String outputFile,String dirPhased, String phasedStart, String phasedEnd,HashMap<String,List<HashMap<String,Object>>> halotypeList){
		System.out.println("---> Generating Haplotype genotypes for " + inputFile);
		
		HashMap<String,HashMap<String,Object>> snpGeno = getPhasedData(homeDir,dirPhased, phasedStart,phasedEnd);

		File file = new File(homeDir + "/" + outputFile + "_"+ inputFile);
		try {
			FileWriter output = new FileWriter(file.getAbsoluteFile());
			BufferedWriter writer = new BufferedWriter(output);
			writer.write("\n");
			for (String hapName : halotypeList.keySet()) {
				HashMap<String,String> genotype = getHaplotypeGenotype(snpGeno,halotypeList.get(hapName));
				if (genotype.get("A1").indexOf("?") != -1 ) {
					writer.write(genotype.get("chromosome") +"\t"+genotype.get("position")+"\t"+ hapName+"\t" +"NA\n");
				}
				else{
					writer.write(genotype.get("chromosome") +"\t"+genotype.get("position")+"\t"+ hapName+"\t" +genotype.get("A1")+"/"+genotype.get("A2")+"\n");
				}
			}
			appendSnpHaplotype(writer,homeDir,inputFile);
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		deleteDirectory(new File(homeDir +"/"+ dirPhased));
	}/*** END processBeagleOutput ***/
	
	/***  appendSnpHaplotype ***/
	private static void appendSnpHaplotype(BufferedWriter writer,String homeDir,String inputFile){
		System.out.println("---> appendSnpHaplotype");
	    try { 
		    /* open file */
			FileReader input = new FileReader(homeDir +"/"+ inputFile);
			BufferedReader fileBuffer = new BufferedReader(input);
			/* read file / append */
			String line = "";
			while (( line = fileBuffer.readLine()) != null) {
				if (!line.equals("")){
					writer.write(line +"\n");
				}
			}
			fileBuffer.close();
	    } 
	    catch (Exception e) {
	    	System.out.println("NO FILE " + homeDir +"/"+  inputFile);
			throw new RuntimeException();
		}
	}/*** END appendSnpHaplotype ***/
	
	/***  getHaplotypeGenotype ***/
	private static HashMap<String,String> getHaplotypeGenotype(HashMap<String,HashMap<String,Object>> snpGeno,List<HashMap<String,Object>> snpList){
		System.out.println("---> getHaplotypeGenotype");
		HashMap<String,String> genotype = new HashMap<String,String>();
		genotype.put("A1", "");
		genotype.put("A2", "");
		genotype.put("chromosome", "");
		genotype.put("position", "");
		for (HashMap<String,Object> snpInfo : snpList) {
			/** QC **/// TODO error/warning log
			if (!genotype.get("chromosome").equals("") && !genotype.get("chromosome").equals(snpInfo.get("chromosome"))){
				System.out.println("PROBLEM CHROMOSOME NUMBER HAPLO "+ snpInfo); 
			}
			else if (genotype.get("chromosome").equals("")){
				genotype.put("chromosome", snpInfo.get("chromosome").toString());
			}
			if (genotype.get("position").equals("") ){
				genotype.put("position", snpInfo.get("position").toString());
			}
			else if ( Integer.parseInt(genotype.get("position")) > Integer.parseInt(snpInfo.get("position").toString()) ){
				System.out.println("PROBLEM POSITION NOT SORTED? "+ snpInfo);
				genotype.put("position", snpInfo.get("position").toString());
			}
			genotype.put("A1", genotype.get("A1") + snpGeno.get(snpInfo.get("dbsnp")).get("A1").toString());
			genotype.put("A2", genotype.get("A2") + snpGeno.get(snpInfo.get("dbsnp")).get("A2").toString());
		}
		return genotype;
	}/*** END getHaplotypeGenotype ***/
	
	/***  getPhasedData ***/
	private static HashMap<String,HashMap<String,Object>> getPhasedData(String homeDir, String dirData, String start, String end ) {
		System.out.println("---> getPhasedData");
		
		HashMap<String,HashMap<String,Object>> hapmapGeno = new HashMap<String,HashMap<String,Object>>();
		File[] listFiles = new File(homeDir + "/"+ dirData).listFiles();
		for (File file : listFiles) {
		    if (file.isFile() && file.getName().startsWith(start) && file.getName().endsWith(end)) {
		        try { 
	        	    /* open file */
		    		FileInputStream input = new FileInputStream(homeDir +"/"+ dirData+"/"+ file.getName());
	        	    GZIPInputStream gzis = new GZIPInputStream(input);
	        	    InputStreamReader xover = new InputStreamReader(gzis);
					BufferedReader fileBuffer = new BufferedReader(xover);
					/* read file / get genotype */
					String line = "";
					while (( line = fileBuffer.readLine()) != null) {
						String[] word = line.split(" ");
						if (!word[0].equals("I")) {
							HashMap<String,Object> genotype = new HashMap<String,Object>();
							genotype.put("A1",(Object) word[2]);
							genotype.put("A2", (Object)  word[3]);
							hapmapGeno.put(word[1],genotype);
						}
					}
					fileBuffer.close();
		        } 
		        catch (Exception e) {
		        	System.out.println("NO FILE " + homeDir +"/"+ dirData+"/"+ file.getName());
					throw new RuntimeException();
				}
		    }
		}
		return hapmapGeno;
	}/*** END getPhasedData ***/

	/*** getHaplotypeInfo 
	 * get a list of haplotype
	 * with list of SNPs order by position
	 * 
	 * ***/
	
	private static void getHaplotypeInfo(HashMap<String,List<HashMap<String,Object>>>halotypeList,HashMap<String,HashMap<Double,Object>> snpPosition,HashMap<String,List<Double>> chromosomePosition){
		System.out.println("---> getHaplotypeInfo  ");

		try {
			/*** Connect to Mongodb TODO ***/
			Mongo mongo = new Mongo("50.112.142.17", 27017);
			DB db = mongo.getDB("genophen30");
			db.authenticate("adm1n", "adm1npw0rd".toCharArray());
			
			/*** Recover hapmap individual names and information ***/
			DBCollection collection = db.getCollection("genetics.rep");
			DBCursor disease = collection.find();
			/** Loop on diseases **/
			while (disease.hasNext()) {
				DBObject dbObject = disease.next();
				@SuppressWarnings (value="unchecked")
				List<Object> snpList = (List<Object>) dbObject.get("snps");
				Iterator<Object> iterator = snpList.iterator();
				/* Loop on SNPs */
				while (iterator.hasNext()) {
					@SuppressWarnings (value="unchecked")
					HashMap<String,Object> doc = (HashMap<String,Object>) iterator.next();
					String snp = doc.get("snp").toString();
					if (snp.split("rs").length == 1 && !halotypeList.containsKey(snp) ){
						List<HashMap<String,Object>> hapSnp = getHaplotypeId(snp);
						halotypeList.put(snp,hapSnp);
						/** add order list of snp to chromosome map **/
						addSnpToList( snpPosition, chromosomePosition,hapSnp);
					}
				}
			}			 
		} catch (UnknownHostException e) {
			e.printStackTrace();
		} catch (MongoException e) {
			e.printStackTrace();
		}
	}/*** END getHaplotypeInfo ***/
	
	
	/*** addSnpToList
	 * add position+ snp info to chr/position hasmap
	 * add Position to chr/position list (the list of position per chromosome will be ordered)
	 * ***/
	private static void  addSnpToList(HashMap<String,HashMap<Double,Object>> snpPosition, HashMap<String,List<Double>> chromosomePosition, List<HashMap<String,Object>> hapSnp ){
		System.out.println("---> addSnpToList  ");
		String chromosome =  hapSnp.get(0).get("chromosome").toString();
		if (!(chromosomePosition.containsKey(chromosome) && snpPosition.containsKey(chromosome))){
			chromosomePosition.put(chromosome,new ArrayList<Double>());
			snpPosition.put(chromosome, new HashMap<Double,Object>());
		}
		for (HashMap<String,Object> snpInfo : hapSnp) {
			Double position =  Double.parseDouble(snpInfo.get("position").toString());
			if ( !snpPosition.get(chromosome).containsKey(position)){
				snpPosition.get(chromosome).put(position,snpInfo);
				chromosomePosition.get(chromosome).add(position);
			}					
		}
		Collections.sort(chromosomePosition.get(chromosome));
	}/*** END addSnpToList ***/
	
	
	/*** getHaplotypeId
	 * return the snp list order by position for haplotype "hap"
	 * ***/
	private static List<HashMap<String,Object>> getHaplotypeId(String hap) {
		System.out.println("---> getHaplotypeId  " + hap);
		
		List<HashMap<String,Object>> snpList = new ArrayList<HashMap<String,Object>>();	
		try {
			/*** Connect to Mongodb TODO ***/
			Mongo mongo = new Mongo("50.112.142.17", 27017);
			DB db = mongo.getDB("genetics");
			db.authenticate("adm1n", "adm1npw0rd".toCharArray());
			
			/*** Recover haplotype information ***/
			DBCollection collection = db.getCollection("idhap");
			BasicDBObject query = new BasicDBObject();
			query.put("dbhap",hap);
			DBCursor hapInfo = collection.find(query).sort( new BasicDBObject( "position" , 1 ));
			/** Loop on snp **/
			while (hapInfo.hasNext()) {
				DBObject dbObject = hapInfo.next();
				@SuppressWarnings (value="unchecked")
				HashMap<String,Object> list = (HashMap<String,Object>) dbObject.toMap();
				snpList.add(list );
			}
		} catch (UnknownHostException e) {
			e.printStackTrace();
		} catch (MongoException e) {
			e.printStackTrace();
		}
	    return snpList;
	}/*** END getHaplotypeId ***/
	
	static public void deleteDirectory(File path) 
	{
	    if (path == null)
	        return;
	    if (path.exists())
	    {
	        for(File f : path.listFiles())
	        {
	            if(f.isDirectory()) 
	            {
	                deleteDirectory(f);
	                f.delete();
	            }
	            else
	            {
	                f.delete();
	            }
	        }
	        path.delete();
	    }
	}
}
