package com;
import java.io.*;

import com.mongodb.*;
import java.util.*; 
import java.lang.reflect.Method;  
import java.net.*;  
import java.util.jar.*;  
import java.util.zip.GZIPInputStream;


public class GzipTest {
	public static void main(String args[]) {
		String path ="/Users/celine/Box Documents/Celine/ScriptInOut/mongodb-update";
		String filename = "PG0000802-BLD.snps.vcf.gz";
		try{
			FileInputStream input = new FileInputStream(path +"/"+filename);
		    GZIPInputStream gzis = new GZIPInputStream(input);
		    InputStreamReader xover = new InputStreamReader(gzis);
			BufferedReader fileBuffer = new BufferedReader(xover);
			/* read file / get genotype */
			String line = "";
			line = fileBuffer.readLine();
			System.out.println(line);
		}
		catch (Exception a){
			
		}
		
    }/*** END MAIN ***/

}