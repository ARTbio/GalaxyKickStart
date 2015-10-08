import java.lang.reflect.Field;

public class PMC
{
	public class ProcessPair
	{
		public Process pr;
		public String command;
	}

	public int cpus;
	public long pid, start;
	public String inDirName;
	public java.io.File finishedFile;
	public java.io.File executingFile;
	public java.io.File waitingFile;

	public java.util.Vector < ProcessPair > executing;
	public java.util.Vector < String > waiting;
	public int waitInterval;
	public Thread sleep;
	
	PMC(int qtd, String in, long jobid)
	{
		pid = 0;
		start = jobid;
		sleep = new Thread();
		sleep.start();
		waitInterval = 1000;
		executing = new  java.util.Vector < ProcessPair > ();
		waiting = new java.util.Vector < String > ();

		inDirName = in;
		cpus = qtd;
		finishedFile = new java.io.File (inDirName + "../" + "/finished.out");
		executingFile = new java.io.File (inDirName + "../" + "/executing.out");
		waitingFile = new java.io.File (inDirName + "../" + "/waiting.out");
		executingFile.delete();
		waitingFile.delete();
	}

	public void append(java.io.File file, String line)
	{
		try
		{
		    java.io.FileWriter fw = new java.io.FileWriter(file.getAbsoluteFile(),true);
		    fw.write(line + "\n");
		    fw.close();
		}
		catch(java.io.IOException ioe)
		{
		    System.err.println("IOException: " + ioe.getMessage());
		}
	}


	public java.util.Vector < String > read(java.io.File file)
	{
		java.util.Vector < String > lines = new java.util.Vector < String > ();
		String strLine;
		
		try
		{
			java.io.FileInputStream fstream = new java.io.FileInputStream(file.getAbsoluteFile());
			java.io.DataInputStream in = new java.io.DataInputStream(fstream);
			java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(in));
			while ((strLine = br.readLine()) != null)
				if (!strLine.trim().isEmpty())
					lines.add(strLine);
			in.close();
		}
		catch (Exception e)
		{
			  System.err.println("Error: " + e.getMessage());
		}
		
		return lines;
	}
	
	public void move(java.io.File file)
	{
	    java.io.File dir = new java.io.File(file.getParent() + "/../");
	    file.renameTo(new java.io.File(dir, file.getName()));
	}
	
	public void update()
	{
		java.io.File inDir = new java.io.File(inDirName);
		java.io.File currentFile, commandFile;
		String [] fileList;
		java.util.Vector < String > lines;
		commandFile = new java.io.File(inDirName + "/command.in");
		
		while (true)
		{
			if (cpus - executing.size() - waiting.size() > 0)
			{
				fileList = inDir.list();
				for (int fileListElement = 0; fileListElement < fileList.length; fileListElement++)
				{
					currentFile = new java.io.File(inDirName + fileList[fileListElement]);
					lines = read(currentFile);
					while(lines.size() > 0 & (cpus - executing.size() > waiting.size()))
					{
						waiting.add(lines.get(0));
						append(waitingFile,lines.get(0));
						lines.remove(0);
					}
					currentFile.delete();
					while(lines.size() > 0)
					{
						append(currentFile,lines.get(0));
						lines.remove(0);
					}
					if (cpus - executing.size() - waiting.size() == 0) break; 
				}
			}

			if (commandFile.exists())
			{
				lines = read(commandFile);
				if (lines.size() > 0)
				{
					execute(lines.get(0));
					System.out.println("command");
				}
				commandFile.delete();
			}
			
			if (checkFinished())
			{
				try
				{
					executingFile.delete();
					executingFile.createNewFile();
					for (int executingElement = 0; executingElement < executing.size(); executingElement++)
						append(executingFile,executing.get(executingElement).command);
				}
				catch (java.io.IOException ioe)
				{
				}
			}
			
			if (submit()) 
			{
				try
				{
					waitingFile.delete();
					waitingFile.createNewFile();
					for (int executingElement = 0; executingElement < waiting.size(); executingElement++)
						append(waitingFile, waiting.get(executingElement));
				}
				catch (java.io.IOException ioe)
				{
				}
			}

			try
			{
				Thread.sleep(waitInterval);
				System.gc();
				executingFile.setLastModified(java.util.Calendar.getInstance().getTimeInMillis());
			}
			catch (InterruptedException ie)
			{
			}
		}
	}
	
	public void execute(String command)
	{
		ProcessBuilder pb;
		
		try
		{
			pb = new ProcessBuilder("bash", "-c", command);
			pb.redirectErrorStream(true);
			pb.environment().put("pmc_job_id", start + "." + pid);
			pb.start();
		}
		catch (java.io.IOException ioe)
		{
		}
	}

	
	public boolean submit()
	{
		boolean result = false;
		ProcessBuilder pb;
		ProcessPair pp;
		Field f;
		
		while (waiting.size() > 0 & executing.size() < cpus)
		{
			pp = new ProcessPair();
			pp.command = waiting.get(0);
			try
			{
				pb = new ProcessBuilder("bash", "-c", pp.command);
				pp.command = start + "." + pid + " - " + pp.command;
				pb.redirectErrorStream(true);
				pb.environment().put("pmc_job_id", start + "." + pid);
				pp.pr = pb.start();

				f = pp.pr.getClass().getDeclaredField("pid");
			    f.setAccessible(true);
			    pp.command = f.getInt(pp.pr) + " - " + pp.command;
				
				executing.add(pp);
				waiting.remove(0);
				append(executingFile,pp.command);
				result = true;
				pid++;
			}
			catch (Exception ioe)
			{
			}
		}

		return result;
	}
	
	public boolean checkFinished()
	{
		boolean freeSlot = false;
		
		for (int slotElement = 0; slotElement < executing.size(); slotElement++ )
		{
			try 
			{
				executing.get(slotElement).pr.exitValue();
				append(finishedFile, executing.get(slotElement).command);
				executing.remove(slotElement);
				freeSlot = true;
			}
			catch (IllegalThreadStateException itse)
			{
			}
		}
		
		return freeSlot;
	}
	
	public static void main (java.lang.String [] arguments)
	{
		PMC pmc = new PMC(Integer.parseInt(arguments[0]), arguments[1], Long.parseLong(arguments[2]));
		pmc.update();
	}
}
