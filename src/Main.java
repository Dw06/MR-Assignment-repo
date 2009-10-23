import javax.swing.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.FlowLayout;

import java.util.*;

import java.io.File;
import java.io.IOException;
//libraries
import jigl.image.*; 
import jigl.image.io.*; 
import jigl.image.ops.Convolve;
import jigl.image.utils.ImageConverter;
import jigl.gui.*;

public class Main implements ActionListener {
	GrayImage original, sobel, hT, maxima,lines;
	JDesktopPane panel;
	JInternalFrame buttons,origFrame, sobelFrame, houghFrame, maxFrame,lineFrame;
	JFrame frame;
	JMenuItem open, save, quit;
	JButton sobelB, houghB, maxB, linesB, detectL;
	
	JImageCanvas hTCan;
	
	boolean finished;
	int thetaMax, rMax, nbrSize;
	int [][] H;
	double [] sinArray, cosArray;

	ArrayList<houghLine> hlines;

	Main() {
		frame = new JFrame("IBV Practicum: The Hough Transform, by David Weterings, 3117480");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(800, 600);
		frame.setJMenuBar(makeFileMenu());
		constructPanels();
		frame.setVisible(true);
	}
	
	public void constructPanels() {
		finished = false;
		original = null; sobel = null; hT = null; maxima = null; lines = null;
		JPanel bt = new JPanel();
		sobelB = new JButton("Sobel Filter");
		houghB = new JButton("Hough Transform");
		maxB = new JButton("Show Maxima");
		linesB = new JButton("Detect lines");
		sobelB.setEnabled(false);
		sobelB.addActionListener(this);
		houghB.setEnabled(false);
		houghB.addActionListener(this);
		maxB.setEnabled(false);
		maxB.addActionListener(this);
		linesB.setEnabled(false);
		linesB.addActionListener(this);
		buttons = new JInternalFrame();
		bt.setLayout(new FlowLayout());
		bt.add(sobelB);
		bt.add(houghB);
		bt.add(maxB);
		bt.add(linesB);
		buttons.add(bt);
        buttons.setSize(165, 170);
        buttons.setVisible(true);
		panel = new JDesktopPane();
		panel.add(buttons);
		frame.add(panel);
		hlines = new ArrayList<houghLine>(); 
	}
	
	public static void main(String[] args) {   
		new Main();
    }
	
	public JMenuBar makeFileMenu() {   
		JMenuBar menubar = new JMenuBar();
        JMenu file;
        
        file = new JMenu("File");
        
        open = new JMenuItem("Open");
        open.getAccessibleContext().setAccessibleDescription("Open");
        open.addActionListener(this);
        file.add(open);
        
        save = new JMenuItem("Save selected");
        save.getAccessibleContext().setAccessibleDescription("Save");
        save.addActionListener(this);
        file.add(save);
        
        quit = new JMenuItem("Quit");
        quit.getAccessibleContext().setAccessibleDescription("Quit");
        quit.addActionListener(this);
        file.add(quit);
        
        menubar.add(file);

        return menubar;
    }

	public void saveImage(Object image, String filename) throws IOException, ImageNotSupportedException	{
		ImageOutputStreamJAI outputJPEG = new ImageOutputStreamJAI(filename + ".jpg");
		outputJPEG.writeJPEG((Image)image);
	}
	
	public void actionPerformed(ActionEvent e) {
			//SHOW SOBEL FILTER
		if(e.getSource() == sobelB) {
			try {
				linesB.setEnabled(false);
				if(!finished) {
					sobel = sobelFilter(original);
				}
				sobelFrame = new JInternalFrame("Sobel Filter");
	   	       	JImageCanvas sobelCanvas = new JImageCanvas(sobel);
	   	       	sobelFrame.add(sobelCanvas);
	   	       	sobelFrame.setSize(sobel.X()+10, sobel.Y()+35);
	   	       	panel.add(sobelFrame);
	   	       	sobelFrame.setLocation(200,50);
	   	       	sobelFrame.setVisible(true);
	   	       	sobelB.setEnabled(false);
	   	       	if(!finished) {
	   	       		houghB.setEnabled(true);
	   	       	}
			} catch(Exception es) {es.printStackTrace();}
			//SHOW HOUGH TRANSFORM
		} else if (e.getSource() == houghB) {
			try {
				if(!finished) {
					hT = houghTransform(sobel);
				}
	   	       	houghFrame = new JInternalFrame("The Hough Transform Image");
	   	       	hTCan = new JImageCanvas(hT);
	   	       	houghFrame.add(hTCan);
	   	       	houghFrame.setSize(hT.X()+10, hT.Y()+35);
	   	       	panel.add(houghFrame);
	   	       	houghFrame.setLocation(200, 80);
	   	       	houghFrame.setVisible(true);
	            houghB.setEnabled(false);
	            
	            if(!finished) {
	            	maxB.setEnabled(true);
	            }
			} catch(Exception es) {es.printStackTrace();}
			//SHOW MAXIMA
		} else if (e.getSource() == maxB) {
			try {
				if(finished) {
					if(houghB.isEnabled()) {
			   	       	houghFrame = new JInternalFrame("The Hough Transform Image");
			   	       	hTCan = new JImageCanvas(hT);
			   	       	houghFrame.add(hTCan);
			   	       	houghFrame.setSize(hT.X()+10, hT.Y()+35);
			   	       	panel.add(houghFrame);
			   	       	houghFrame.setLocation(200, 80);
			   	       	houghFrame.setVisible(true);
			            houghB.setEnabled(false);
					} 
				} else {
					maxima = findMaxima(hT);
				}
				houghFrame.setTitle("Hough + maxima");
				hTCan.setImage(maxima); 
				hTCan.repaint();
	   	       	maxB.setEnabled(false);
	   	       	if(!finished) {
	   	       		linesB.setEnabled(true);
	   	       	}
			} catch(Exception es) {es.printStackTrace();}
			//SHOW DETECTED LINES
		} else if (e.getSource() == linesB) {
			try {
				linesB.setEnabled(false);
				if(!maxB.isEnabled() && !finished) {
					sobelB.setEnabled(true);
					houghB.setEnabled(true);
					maxB.setEnabled(true);
					sobel = sobelFilter(original);
					hT = houghTransform(sobel);
					maxima = findMaxima(hT);
				}
				GrayImage lines = visualizeLines((GrayImage)original.copy());
	            lineFrame = new JInternalFrame("Detected Lines");
	            JImageCanvas lineCan = new JImageCanvas(lines);
	            lineFrame.add(lineCan);
	            lineFrame.setSize(lines.X()+10, lines.Y()+35);
	            panel.add(lineFrame);
	            lineFrame.setLocation(200, 110);
	            lineFrame.setVisible(true);
	            finished = true;
			} catch(Exception es) {es.printStackTrace();}
			//OPEN
		} else if (e.getSource() == open) {
			File file;
			JFileChooser openFile = new JFileChooser();
        	openFile.showOpenDialog(new JFrame());
        	file = openFile.getSelectedFile();
        	if(file != null) {
        		try {
        			//remove and rebuild everything
        			panel.removeAll();
        			constructPanels();
        			ImageInputStreamJAI i = new ImageInputStreamJAI(file.getPath());
        			Image image = i.read();
        			//zet het plaatje om naar GrayImage
        			image = ImageConverter.toGray(image);
        			original = (GrayImage) image;
        			origFrame = new JInternalFrame("Original Image");
        			JImageCanvas oriCan = new JImageCanvas(original);
        			origFrame.add(oriCan);
        			origFrame.setSize(original.X()+10, original.Y()+35);
        			origFrame.setLocation(200,20);
        			panel.add(origFrame);
        			origFrame.setVisible(true);
        			sobelB.setEnabled(true);
        			linesB.setEnabled(true);
        		 }
        		 catch(Exception ex){
        			 ex.printStackTrace();
        		 }
        	}
        	//SAVE
		} else if(e.getSource() == save) {
			if(origFrame != null) {
				if(origFrame.isSelected()) {
					try {
						saveImage(original, "original");
					} catch (Exception e1) {
						e1.printStackTrace();
					}
				}
			} 
			
			if(sobelFrame != null) {
				if(sobelFrame.isSelected()) {
					try {
						saveImage(sobel, "sobel");
					} catch (Exception e1) {
						e1.printStackTrace();
					}
				}
			} 
			
			if(houghFrame != null) {
				if(houghFrame.isSelected()) {
					try {
						saveImage(hT, "hough");
					} catch (Exception e1) {
						e1.printStackTrace();
					}
				}
			} 
			if(lineFrame != null) {
				if(lineFrame.isSelected()) {
					try {
						saveImage(lines, "detectedlines");
					} catch (Exception e1) {
						e1.printStackTrace();
					}
				}
			}
			//QUIT
		} else if (e.getSource() == quit) {
			System.exit(0);
		}
		
	}
	
	public void precalculations() {	
		for(double i=-90; i<thetaMax-90;i++) {
			double rad = Math.toRadians(i);
			sinArray[(int) (i+90)] = Math.sin(rad);
			cosArray[(int) (i+90)] = Math.cos(rad);
		}
	}
	
	//doe hough transformatie
	public GrayImage houghTransform(GrayImage img) {
		thetaMax = 270;
		rMax = (int) Math.sqrt(img.X()*img.X() + img.Y()*img.Y());

		sinArray = new double[thetaMax];
		cosArray = new double[thetaMax];
		precalculations();
		H = new int[thetaMax][rMax];
		
		int maxValue = Integer.MIN_VALUE;
		int minValue = Integer.MAX_VALUE;
		int pixelValue = 0;
		int r = 0;
		//voor alle object pixels (pixelValue > 0)
		//voor iedere theta, bereken r. Als r geldig is voeg hem toe aan de Hough array
		for(int x=0; x<img.X();x++)
			for(int y=0; y<img.Y();y++)
				for(int ai=-90; ai<thetaMax-90; ai++)	{	
					pixelValue = img.get(x,y);
					if(pixelValue > 0) {
    					r = (int)(x * cosArray[ai+90] + y*sinArray[ai+90]);

						if(r>0 && r<rMax) {
							pixelValue = H[(int) (ai+90)][r] + pixelValue;
							H[(int) (ai+90)][r] = pixelValue;
							if(maxValue < pixelValue) {
								maxValue = pixelValue;
							}
							if(minValue > pixelValue) {
								minValue = pixelValue;
							}
						}
					}
				}
		
		GrayImage result = new GrayImage(thetaMax, rMax);
		//rescaleImage(result, maxValue, minValue);
		result.byteSize();
		return result;	
	}
	
	//rescale een image terug naar range 0..255
	public void rescaleImage(GrayImage gi, int max, int min) {
		int range = max - min;
		int deling = range/255;
		for(int x=0; x<gi.X(); x++) {
			for(int y=0; y<gi.Y(); y++) {
				gi.set(x, y, (H[x][y] - min)/deling);
			}
		}
	}
	
	public GrayImage sobelFilter(Image grayImage) throws ImageNotSupportedException,InvalidKernelException {
		GrayImage temp = (GrayImage) grayImage.copy();
		Convolve cv1 = new Convolve(new ImageKernel(jigl.image.ImageKernel.SOBEL_X));
		Convolve cv2 = new Convolve(new ImageKernel(jigl.image.ImageKernel.SOBEL_Y));
		GrayImage xD = (GrayImage) cv1.apply(grayImage);
		GrayImage yD = (GrayImage) cv2.apply(grayImage);

		int width, height, value,xValue, yValue;
		width = grayImage.X();
		height = grayImage.Y();
		
		for(int x = 0; x<width; x++) {
			for(int y=0; y<height; y++) {
				xValue = xD.get(x,y); 
				yValue = yD.get(x,y);
		    	value = (int) Math.sqrt( (xValue*xValue) + (yValue*yValue));
		    	temp.set(x,y,value);
		    }
		}
		//rescale
		temp.byteSize();
		return temp; 
	}

	public GrayImage findMaxima(GrayImage hough) {
		GrayImage temp = (GrayImage) hough.copy();
		
		nbrSize = 4;
		int max_value = Integer.MIN_VALUE;
		//vind de maximale waarde
		for(int t = 0; t < hough.X(); t++)
			for(int r = 0; r < hough.Y(); r++)
				if (H[t][r] > max_value)
					max_value = H[t][r];
			
		boolean maxima;
		//check de locale maxima, als je er een hebt gevonden voeg hem toe.
		for (int t = 4; t < hough.X() - 5; t++) {
			for (int r = 4; r < hough.Y() - 5; r++)	{
				if (H[t][r] > 0.3*max_value) {
					maxima = true;
					for(int dr = r-4; dr < r+5; dr++) {
						for(int dt = t-5; dt < t+5; dt++) {
							if (H[dt][dr] > H[t][r]) {
								maxima = false;
							}
						}
					}
					if(maxima) {
						hlines.add(new houghLine(t,r,H[t][r]));
					}
				}
			}
		} 
		visualizeMaxima(temp);
		temp.byteSize();
		return temp;
	}
	
	//show Maxima
	public void visualizeMaxima(GrayImage img) {
		Iterator<houghLine> i = hlines.iterator();
		while(i.hasNext()) {
			houghLine hL = i.next();
			hL.drawMax(img);
		}
	}
	
	//show detected lines
	public GrayImage visualizeLines(GrayImage img) {
		Iterator<houghLine> i = hlines.iterator();
		while(i.hasNext()) {
			houghLine hL = i.next();
			hL.draw(img,rMax);
		}
		return img;
	}
}
