import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
import sys





def read_input():
	#read input
	nr_points=int(input("Enter number of points:")) 
	nr_meas=int(input("Enter number of measurements:"))

	#array with 4 rows: 
	#0 point from where measurements goes
	#1 point to where measurement goes
	#2 value of measurement
	#3 relative variance 

	arr_in=np.zeros([nr_meas,4], dtype=float)
	i=0;
	while i<len(arr_in):
		try:
			line_in=raw_input("Enter measurement (4 Elements per line):").strip().split()
			for j in range(0, 4):
				arr_in[i][j]=line_in[j]
			i=i+1;
			
		except IndexError:
			print("Oops!  That was not a valid input.  Try again...")

		
	return {'nr_meas': nr_meas, 'nr_points':nr_points, 'arr_in':arr_in}




def calc(nr_points, nr_meas, arr_in):

	print("The input: " + str(arr_in))

	L=np.zeros([nr_meas,1], dtype=float)
	index=0
	for iii in arr_in:
		L[index]=iii[2]
		index=index+1

	#least square adjustment for getting real values from measurement
	#point 1 is leveled to 0

	#creating A matrix
	arr_a=np.zeros([nr_meas,nr_points-1], dtype=int)

	index=0
	for ii in arr_in:
		start=ii[0]
		end=ii[1]
		if start!=1:
			s=start-2
			arr_a[index][s]=-1  #minus 2 because point 1 is omitted and array starts with 0 not 1
		if end!=1:
			e=end-2
			arr_a[index][e]=1
		index=index+1

	arr_aT=np.transpose(arr_a) # creating A^T 
	aT_a=arr_aT.dot( arr_a) # A^T*a

	aT_L=arr_aT.dot(L) # A^T * L

	x_dach=np.linalg.solve(aT_a, aT_L) # A^T * A* X' =A^T*L , solving the system of linear scalar equations
	print("The solution for the real hights of point 2 to " +str(nr_points)+" in reverence to point 1: \n" + str(x_dach))

	L_dach= arr_a.dot(x_dach) # A * X'
	print("The corrected measurements are: \n"+str(L_dach))



	#calculation of standard deviation
	#P is the identity matrix here, therefore omitted

	Q=np.zeros([nr_meas, nr_meas], dtype=float)

	for jj in range(0, len(arr_in)):
		Q[jj][jj]=arr_in[jj][3]
	#print(Q)

	p=np.linalg.inv(Q)
	#print(p)

	v=L_dach-L
	vT=np.transpose(v)
	if(nr_meas-(nr_points-1)==0):
		quit()

	sigma2=float(vT.dot(p).dot(v)/(nr_meas-(nr_points-1))) #determinante stattdessen?
	sigma=np.sqrt(abs(sigma2))
	print("The standart deviation sigma is: "+ str(sigma))


	Exx=sigma2*np.linalg.inv(arr_aT.dot(p).dot(arr_a)) #precision for heights. Exx= sigma^2 *( A^T * P * A)^(-1)
	out=""
	for i in range(0, len(Exx)):
		out=out+"sigma for point "+ str(i+2)+": "+str(np.sqrt(abs(Exx[i][i])))+"\n"
	print(out)

	Ell=arr_a.dot(Exx).dot(arr_aT) #precisicion for measurments. Ell=A * Exx * A^T
	out=""
	for i in range(0, len(Ell)):
		out=out+"sigma for meansurement from "
		out=out+str( int(arr_in[i][0]))+" to "+str(int(arr_in[i][1]))
		out=out+": "+str(np.sqrt(abs(Ell[i][i])))+"\n"
	print(out)


def gaussian(arr_in):
	for i in arr_in:
	 	i[2]=np.random.normal(loc=i[2], scale=0.2, size=None)
	return arr_in

def main():
	
	res=read_input()
	nr_meas=res['nr_meas']
	nr_points=res['nr_points']
	arr_in=res['arr_in']
		
	calc(nr_points, nr_meas, arr_in)

	
main()
