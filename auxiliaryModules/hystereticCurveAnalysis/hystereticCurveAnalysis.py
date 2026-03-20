#-*-coding: UTF-8-*-
#####Units: Length-m, Force-kN, mass-ton, Stress-kpa(10e-3MPa), g=9.81m/s2
#####Units: Length-mm, Force-N, mass-ton, Stress-Mpa, g=9810mm/s2 pho=ton/mm3
########################################################################################################################
#  Author: Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo
#  E-mail: guojj@tongji.edu.cn/guojj_ce@163.com
#  Environemet: Successfully executed in python 3.11
#  Date: 2024-03-22
########################################################################################################################
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import openseespy.opensees as ops
########################################################################################################################
########################################################################################################################
class HystereticCurveAnalysis():
	"""
	--------------------------------------------------------------------------------------------------------------------
	A class for hysteretic curve analyses (version:0.6.0)
	Environemet: Successfully executed in python 3.11
	Date: 2024-03-22
	--------------------------------------------------------------------------------------------------------------------
	    ** ********************************************************************************  **
	    ** (C) Copyright 2024, School of Civil Engineering, Beijing Jiaotong University      **
	    ** All Rights Reserved.                                                              **
	    **                                                                                   **
	    ** Commercial use of this program is strictly prohibited.                            **
	    **                                                                                   **
	    ** Developed by:                                                                     **
	    **   Junjun Guo,Beijing Jiaotong University. https://github.com/Junjun1guo           **
	    **   jjguo2@bjtu.edu.cn/guojj_ce@163.com                                             **
	    ** ********************************************************************************  **
	"""
	def __init__(self,xDataList,yDataList):
		"""
        ----------------------------------------------------------------------------------------------------------------
        Initialize the class
        ------------------------------------------
        Inputs:
        xDataList,yDataList(list)-The x and y values of the hysteretic curve data
        ----------------------------------------------------------------------------------------------------------------
		"""
		self.xValueList,self.yValueList=xDataList,yDataList
		self.loopInfoDict={} ###---loopNum:[startIndex,endIndex]
		self.tensionOnly=True

	def plotHystereticCurve(self,saveFig=False,xlabel="x",ylabel='y',title='Hysteretic curve',multiColors=False):
		"""
        Plot the original hysteratic curve
        ------------------------------------------
        Inputs:
        saveFig(bool)-Save the plot to eps figure, default not saved
        xDataList,yDataList(list)-The x and y values of the hysteretic curve data
        multiColors(bool)-Whether distinguish different loops with different colors, default value is False
		"""
		ax=plt.subplot(1,1,1)
		if not multiColors:
			plt.plot(self.xValueList,self.yValueList, ls='-', color='b', lw=1.5)
		else:
			self._hystereticLoopDecomp()
			legendList=list(self.loopInfoDict.keys())
			legendList=[str(each) for each in legendList]
			for key,value in self.loopInfoDict.items():
				start,end=value
				plt.plot(self.xValueList[start:end],self.yValueList[start:end],ls='-',lw=1.5)
			ax.legend(labels=legendList)
		plt.xlabel(f'{xlabel}', fontsize=15)
		plt.ylabel(f'{ylabel}', fontsize=15)
		plt.title(f'{title}', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		if saveFig:
			plt.savefig(title+".eps")
			plt.savefig(title + ".jpg")
		plt.show()

	def skeletonCurve(self,saveData=False,saveFig=False):
		"""
		Get and plot the skeleton curve of the hysteretic data
        ------------------------------------------
        Inputs:
        saveData(bool)-Whether save the skeleton curve data, default value is False
        saveFig(bool)-Whether save the skeleton curve figure, default value is False
		"""
		posIndex,negIndex=[],[]
		self._hystereticLoopDecomp()
		for eachI in self.loopInfoDict.values():
			posValue=max(self.yValueList[eachI[0]:eachI[1]])
			indexValue=[index for index,value in enumerate(self.yValueList[eachI[0]:eachI[1]])
						if value==posValue][0]+eachI[0]
			posIndex.append(indexValue)
		if not self.tensionOnly:
			for eachI in self.loopInfoDict.values():
				negValue = min(self.yValueList[eachI[0]:eachI[1]])
				indexValue = [index for index, value in enumerate(self.yValueList[eachI[0]:eachI[1]])
							  if value == posValue][0] + eachI[0]
				negIndex.append(indexValue)
		if len(negIndex)>0:
			negIndex.reverse()
		plt.subplot(1, 1, 1)
		plt.plot(self.xValueList, self.yValueList, ls='-', color='b', lw=1)
		plt.xlabel(f'x', fontsize=15)
		plt.ylabel(f'y', fontsize=15)
		plt.title(f'Hysteretic curve', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		negX, negY = [], []
		if len(negIndex)>0:
			[[negX.append(self.xValueList[each]),negY.append(self.yValueList[each])] for each in negIndex]
			for item in range(len(self.xValueList)):
				indexI=negIndex[-1]+item
				if self.yValueList[indexI]>0:
					break
				else:
					negX.append(self.xValueList[indexI])
					negY.append(self.yValueList[indexI])
			plt.plot(negX, negY, ls='-', color='r', lw=1.5)
		posX,posY=[],[]
		for i1 in range(posIndex[0]):
			posX.append(self.xValueList[i1])
			posY.append(self.yValueList[i1])
		for i2 in range(len(posIndex)-1):
			posX.append(self.xValueList[posIndex[i2+1]])
			posY.append(self.yValueList[posIndex[i2+1]])
		plt.plot(posX, posY, ls='-', color='r', lw=1.5)

		if saveData:
			xData=negX+posX
			yData=negY+posY
			saveDataValue=[]
			for eachx,eachy in zip(xData,yData):
				saveDataValue.append([eachx,eachy])
			np.savetxt("skeletonData.txt",saveDataValue, fmt="%.6f %.6f")
		if saveFig:
			plt.savefig("skeleton.eps")
			plt.savefig("skeleton.jpg")
		plt.show()

	def skeletonCurve_doubleDirection(self):
		"""
		Get and plot the skeleton curve of double direction hysteretic data
        ------------------------------------------
		"""
		min_points=50
		x, y =self.xValueList,self.yValueList
		print(f"Total points: {len(x)}, Disp range: [{x.min():.2f}, {x.max():.2f}]")
		zero_tol = 0.02 * np.max(np.abs(x))
		# ---- 2. Find zero-crossing regions ----
		near_zero = np.abs(x) <= zero_tol
		diff_nz = np.diff(near_zero.astype(int))
		enters = np.where(diff_nz == 1)[0] + 1
		leaves = np.where(diff_nz == -1)[0] + 1
		if near_zero[0]:
			enters = np.concatenate([[0], enters])
		if near_zero[-1]:
			leaves = np.concatenate([leaves, [len(x)]])
		n = min(len(enters), len(leaves))
		zero_pts = sorted(set((enters[:n] + leaves[:n]) // 2))
		if not zero_pts or zero_pts[0] > min_points:
			zero_pts.insert(0, 0)
		# ---- 3. Extract peak points from each half-cycle ----
		pos_peaks_disp = []
		pos_peaks_force = []
		neg_peaks_disp = []
		neg_peaks_force = []
		for i in range(len(zero_pts) - 1):
			s, e = zero_pts[i], zero_pts[i + 1]
			if e - s < min_points:
				continue
			xs = x[s:e + 1]
			ys = y[s:e + 1]
			# Find the index of the peak absolute displacement in this half-cycle
			idx_peak = np.argmax(np.abs(xs))
			peak_disp = xs[idx_peak]
			peak_force = ys[idx_peak]
			if peak_disp > 0:
				pos_peaks_disp.append(peak_disp)
				pos_peaks_force.append(peak_force)
			else:
				neg_peaks_disp.append(peak_disp)
				neg_peaks_force.append(peak_force)
		# ---- 4. Build skeleton curve arrays ----
		# Add origin point (0, 0) to both branches for a complete skeleton
		# Positive branch: sort by displacement ascending
		pos_peaks_disp = np.array(pos_peaks_disp)
		pos_peaks_force = np.array(pos_peaks_force)
		neg_peaks_disp = np.array(neg_peaks_disp)
		neg_peaks_force = np.array(neg_peaks_force)
		# Sort positive branch by displacement (ascending)
		if len(pos_peaks_disp) > 0:
			sort_pos = np.argsort(pos_peaks_disp)
			pos_peaks_disp = pos_peaks_disp[sort_pos]
			pos_peaks_force = pos_peaks_force[sort_pos]
			# Remove duplicate displacement levels: keep the one with max force
			unique_pos_disp = []
			unique_pos_force = []
			visited = set()
			# Round to a tolerance to group similar displacement amplitudes
			disp_round = np.round(pos_peaks_disp, decimals=0)
			for d_val in np.unique(disp_round):
				mask = disp_round == d_val
				idx_max_f = np.argmax(pos_peaks_force[mask])
				unique_pos_disp.append(pos_peaks_disp[mask][idx_max_f])
				unique_pos_force.append(pos_peaks_force[mask][idx_max_f])
			pos_peaks_disp = np.array(unique_pos_disp)
			pos_peaks_force = np.array(unique_pos_force)
		# Sort negative branch by displacement (descending, i.e., from 0 toward most negative)
		if len(neg_peaks_disp) > 0:
			sort_neg = np.argsort(neg_peaks_disp)[::-1]
			neg_peaks_disp = neg_peaks_disp[sort_neg]
			neg_peaks_force = neg_peaks_force[sort_neg]
			# Remove duplicate displacement levels: keep the one with min (most negative) force
			unique_neg_disp = []
			unique_neg_force = []
			disp_round = np.round(neg_peaks_disp, decimals=0)
			for d_val in np.unique(disp_round)[::-1]:
				mask = disp_round == d_val
				idx_min_f = np.argmin(neg_peaks_force[mask])
				unique_neg_disp.append(neg_peaks_disp[mask][idx_min_f])
				unique_neg_force.append(neg_peaks_force[mask][idx_min_f])
			neg_peaks_disp = np.array(unique_neg_disp)
			neg_peaks_force = np.array(unique_neg_force)
		# Prepend origin to both branches
		skeleton_pos = np.column_stack([
			np.concatenate([[0.0], pos_peaks_disp]),
			np.concatenate([[0.0], pos_peaks_force])
		])
		skeleton_neg = np.column_stack([
			np.concatenate([[0.0], neg_peaks_disp]),
			np.concatenate([[0.0], neg_peaks_force])
		])
		print(f"Skeleton curve: {len(skeleton_pos)} positive points, {len(skeleton_neg)} negative points")
		# ---- 5. Plot ----
		fig, ax = plt.subplots(figsize=(12, 7))
		# Raw hysteresis in gray
		ax.plot(x, y, color='gray', linewidth=0.3, alpha=0.5, label='Hysteresis Curve')
		# Skeleton curve in blue
		ax.plot(skeleton_pos[:, 0], skeleton_pos[:, 1],
				'b-o', linewidth=2.0, markersize=6, markerfacecolor='blue',
				markeredgecolor='white', markeredgewidth=0.8, zorder=5,
				label='Skeleton Curve (Positive)')
		ax.plot(skeleton_neg[:, 0], skeleton_neg[:, 1],
				'b-s', linewidth=2.0, markersize=6, markerfacecolor='blue',
				markeredgecolor='white', markeredgewidth=0.8, zorder=5,
				label='Skeleton Curve (Negative)')
		# Mark peak points
		ax.plot(skeleton_pos[1:, 0], skeleton_pos[1:, 1], 'ro', ms=4, zorder=6, alpha=0.6)
		ax.plot(skeleton_neg[1:, 0], skeleton_neg[1:, 1], 'ro', ms=4, zorder=6, alpha=0.6)
		# Annotate peak values
		for i in range(1, len(skeleton_pos)):
			ax.annotate(f'({skeleton_pos[i, 0]:.1f}, {skeleton_pos[i, 1]:.1f})',
						xy=(skeleton_pos[i, 0], skeleton_pos[i, 1]),
						xytext=(5, 8), textcoords='offset points',
						fontsize=7, color='blue', alpha=0.8)
		for i in range(1, len(skeleton_neg)):
			ax.annotate(f'({skeleton_neg[i, 0]:.1f}, {skeleton_neg[i, 1]:.1f})',
						xy=(skeleton_neg[i, 0], skeleton_neg[i, 1]),
						xytext=(5, -12), textcoords='offset points',
						fontsize=7, color='blue', alpha=0.8)

		ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
		ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
		ax.set_xlabel('Displacement', fontsize=12)
		ax.set_ylabel('Force', fontsize=12)
		ax.set_title('Hysteresis Curve with Skeleton (Backbone) Curve', fontsize=14)
		ax.legend(loc='upper left', fontsize=10)
		ax.grid(True, alpha=0.3)
		plt.tight_layout()
		plt.show()
		return skeleton_pos, skeleton_neg

	def _hystereticLoopDecomp(self):
		"""
		Decompositon of the hysteretic curves
        ------------------------------------------
        Inputs:
		"""
		loopInfoList=[]
		startIndex=0
		xValue,yValue=self.xValueList,self.yValueList
		xValue[0]=0.0
		if min(xValue)<-0.01*max(xValue): ###---for tensile and compresive loops
			self.tensionOnly=False
			for i1 in range(startIndex, len(xValue)-1):
				if ((xValue[i1+1]-xValue[i1])>0)and(yValue[i1]*yValue[i1+1]<=0):
					loopInfoList.append(i1)
		else:###---for tension only loops
			for i1 in range(1,len(xValue)-1):
				if (xValue[i1]-xValue[i1-1]<0)and(xValue[i1+1]-xValue[i1]>0):
					loopInfoList.append(i1)
		loopInfoList=[0]+loopInfoList+[len(xValue)-1]
		for num in range(len(loopInfoList)-1):
			self.loopInfoDict[num+1]=[loopInfoList[num],loopInfoList[num+1]]

	def plotLoop_singleDirection(self,loopNumber,saveData=False,saveFig=False,dottedLine=True):
		"""
		Plot each hysteretic loop
        ------------------------------------------
        Inputs:
        loopNumber(int)-The hysteretic loop number
        saveData(bool)-Whether save the hysteretic loop data, default value is False
        saveFig(bool)-Whether save the hysteretic loop figure, default value is False
		"""
		self._hystereticLoopDecomp()
		start,end=self.loopInfoDict[loopNumber]
		plt.subplot(1, 1, 1)
		if dottedLine:
			plt.plot(self.xValueList[start:end], self.yValueList[start:end], '.', color='b',markersize=5)
		else:
			plt.plot(self.xValueList[start:end], self.yValueList[start:end], ls='-', color='b', lw=1)
		plt.xlabel(f'x', fontsize=15)
		plt.ylabel(f'y', fontsize=15)
		plt.title(f'loop-{loopNumber}', fontsize=15)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.grid(True)
		if saveFig:
			plt.savefig(f"loop-{loopNumber}.eps")
			plt.savefig(f"loop-{loopNumber}.jpg")
		if saveData:
			saveDataValue=[]
			for eachx,eachy in zip(self.xValueList[start:end], self.yValueList[start:end]):
				saveDataValue.append([eachx,eachy])
			np.savetxt(f"loop-{loopNumber}.txt",saveDataValue, fmt="%.6f %.6f")
		plt.show()

	def plotLoop_doubleDirection(self,loopNumber,saveData,saveFig,dottedLine,zero_tol,min_points):
		"""
		Plot each hysteretic loop
        ------------------------------------------
        Inputs:
        loopNumber(int)-The hysteretic loop number
        saveData(bool)-Whether save the hysteretic loop data, default value is False
        saveFig(bool)-Whether save the hysteretic loop figure, default value is False
        zero_tol(float)-The zero tolerance, default value is None
        min_points(int)-minimum data points per half-cycle
		"""
		x_raw, y_raw =self.xValueList,self.yValueList
		print(f"Total points: {len(x_raw)}, Disp range: [{x_raw.min():.2f}, {x_raw.max():.2f}]")
		if zero_tol is None:
			zero_tol = 0.02 * np.max(np.abs(x_raw))
		# ---- 2. Find zero-crossing regions ----
		near_zero = np.abs(x_raw) <= zero_tol
		diff_nz = np.diff(near_zero.astype(int))
		enters = np.where(diff_nz == 1)[0] + 1
		leaves = np.where(diff_nz == -1)[0] + 1
		if near_zero[0]:
			enters = np.concatenate([[0], enters])
		if near_zero[-1]:
			leaves = np.concatenate([leaves, [len(x_raw)]])
		n = min(len(enters), len(leaves))
		zero_pts = sorted(set((enters[:n] + leaves[:n]) // 2))
		if not zero_pts or zero_pts[0] > min_points:
			zero_pts.insert(0, 0)
		# ---- 3. Extract half-cycles ----
		half_cycles = []
		for i in range(len(zero_pts) - 1):
			s, e = zero_pts[i], zero_pts[i + 1]
			if e - s < min_points:
				continue
			xs, ys = x_raw[s:e + 1], y_raw[s:e + 1]
			peak = xs[np.argmax(np.abs(xs))]
			half_cycles.append({'x': xs.copy(),'y': ys.copy(),'peak_disp': peak,'direction': 'positive' if peak > 0
			else 'negative'})
		# ---- 4. Pair into full loops and calculate damping ratio ----
		full_loops = []
		i = 0
		while i < len(half_cycles) - 1:
			h1, h2 = half_cycles[i], half_cycles[i + 1]
			if h1['direction'] != h2['direction']:
				x_full = np.concatenate([h1['x'], h2['x'][1:]])
				y_full = np.concatenate([h1['y'], h2['y'][1:]])
				# Get positive and negative peak displacement (max/min x in the loop)
				# and their corresponding forces.
				idx_pos_peak_in_loop = np.argmax(x_full)
				idx_neg_peak_in_loop = np.argmin(x_full)
				pos_peak_disp = x_full[idx_pos_peak_in_loop]
				pos_peak_force = y_full[idx_pos_peak_in_loop]
				neg_peak_disp = x_full[idx_neg_peak_in_loop]
				neg_peak_force = y_full[idx_neg_peak_in_loop]
				# Energy dissipated (area of the loop)
				energy_dissipated = abs(np.trapz(y_full, x_full))
				# Equivalent elastic strain energy (Es)
				# Es = 0.5 * (pos_peak_disp * pos_peak_force + |neg_peak_disp * neg_peak_force|)
				# This represents the area of the two triangles formed by the origin and the peak points.
				elastic_strain_energy_es = 0.5 * (pos_peak_disp * pos_peak_force + abs(neg_peak_disp * neg_peak_force))
				# Calculate equivalent damping ratio
				eq_damping_ratio = np.nan  # Initialize as NaN for cases where calculation is not possible
				if elastic_strain_energy_es > 1e-12:  # Avoid division by zero or very small numbers
					eq_damping_ratio = energy_dissipated / (2 * np.pi * elastic_strain_energy_es)
				full_loops.append({
					'id': len(full_loops) + 1,'x': x_full,'y': y_full,'pos_peak_disp': pos_peak_disp,
					'pos_peak_force': pos_peak_force,'neg_peak_disp': neg_peak_disp,'neg_peak_force': neg_peak_force,
					'energy_dissipated': energy_dissipated,'elastic_strain_energy_es': elastic_strain_energy_es,
					'eq_damping_ratio': eq_damping_ratio,'num_points': len(x_full)})
				i += 2
			else:
				i += 1

		print(f"Extracted {len(full_loops)} complete hysteresis loops\n")

		# ---- 5. Plot ----
		# --- Plot a specific loop ---
		if loopNumber< 1 or loopNumber > len(full_loops):
			print(f"Invalid loop_id. Please enter 1 ~ {len(full_loops)}")
			return full_loops
		lp = full_loops[loopNumber - 1]
		fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # Increased figure size for more info
		# Left: highlight in full data
		axes[0].plot(x_raw, y_raw, 'gray', linewidth=0.3, alpha=0.3, label='Raw data')
		axes[0].plot(lp['x'], lp['y'], 'r-', linewidth=1.5, label=f'Loop #{loopNumber}')
		axes[0].plot(lp['x'][0], lp['y'][0], 'go', ms=8, zorder=5, label='Start')
		axes[0].plot(lp['x'][-1], lp['y'][-1], 'bs', ms=8, zorder=5, label='End')
		axes[0].set_xlabel('Displacement')
		axes[0].set_ylabel('Force')
		axes[0].set_title(f'Loop #{loopNumber} in Full Curve')
		axes[0].legend()
		axes[0].grid(True, alpha=0.3)
		# Right: detailed view
		axes[1].plot(lp['x'], lp['y'], 'b-', linewidth=2)
		axes[1].fill(lp['x'], lp['y'], color='blue', alpha=0.06)
		# Plot peak points
		pk_pos_idx = np.argmax(lp['x'])
		pk_neg_idx = np.argmin(lp['x'])
		axes[1].plot(lp['x'][pk_pos_idx], lp['y'][pk_pos_idx], 'r^', ms=10, zorder=5,
					 label=f"+Peak (D,F): ({lp['pos_peak_disp']:.2f}, {lp['pos_peak_force']:.2f})")
		axes[1].plot(lp['x'][pk_neg_idx], lp['y'][pk_neg_idx], 'bv', ms=10, zorder=5,
					 label=f"-Peak (D,F): ({lp['neg_peak_disp']:.2f}, {lp['neg_peak_force']:.2f})")
		# Optional: Plot triangles for visualizing elastic strain energy
		axes[1].plot([0, lp['pos_peak_disp'], lp['pos_peak_disp'], 0],
					 [0, 0, lp['pos_peak_force'], lp['pos_peak_force']], 'r:', alpha=0.5, linewidth=1)
		axes[1].plot([0, lp['neg_peak_disp'], lp['neg_peak_disp'], 0],
					 [0, 0, lp['neg_peak_force'], lp['neg_peak_force']], 'g:', alpha=0.5, linewidth=1)
		info = (f"+Peak Disp: {lp['pos_peak_disp']:.4f}\n"
				f"-Peak Disp: {lp['neg_peak_disp']:.4f}\n"
				f"+Peak Force: {lp['pos_peak_force']:.4f}\n"
				f"-Peak Force: {lp['neg_peak_force']:.4f}\n"
				f"Energy Dissipated (Ed): {lp['energy_dissipated']:.4f}\n"
				f"Elastic Strain Energy (Es): {lp['elastic_strain_energy_es']:.4f}\n"
				f"Equivalent Damping Ratio (eq): {lp['eq_damping_ratio']:.4f}\n"
				f"Points: {lp['num_points']}")
		axes[1].text(0.05, 0.95, info, transform=axes[1].transAxes, fontsize=10,
					 va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
		axes[1].set_xlabel('Displacement')
		axes[1].set_ylabel('Force')
		axes[1].set_title(f'Loop #{loopNumber} Detail')
		axes[1].legend(fontsize='small')
		axes[1].grid(True, alpha=0.3)
		axes[1].axhline(0, color='gray', linewidth=0.5, linestyle='--')
		axes[1].axvline(0, color='gray', linewidth=0.5, linestyle='--')
		plt.tight_layout()
		plt.show()
		# ---- 6. Summary ----
		print(f"\n{'No.':<4} {'+Disp':<10} {'+Force':<10} {'-Disp':<10} {'-Force':<10} {'Ed':<12} {'Es':<12} {'Xi_eq':<12} {'Points':<8}")
		print("-" * 90)
		for lp in full_loops:
			print(f"{lp['id']:<4} {lp['pos_peak_disp']:<10.3f} {lp['pos_peak_force']:<10.3f} "
				  f"{lp['neg_peak_disp']:<10.3f} {lp['neg_peak_force']:<10.3f} "
				  f"{lp['energy_dissipated']:<12.3f} {lp['elastic_strain_energy_es']:<12.3f} "
				  f"{lp['eq_damping_ratio']:<12.4f} {lp['num_points']:<8}")
		return full_loops

	def yAxisDataTranslation(self,startX:float,endX:float):
		"""
		By adjusting the y-direction data of the initial curve to make the data symmetric within the starting interval.
        ----------------------------------------------------------------------------------------------------------------
        Inputs:
        		startX(float)-The start value of the range
        		endX(float)-The end value of the range
		"""
		mask = (self.xValueList >=startX) & (self.xValueList <=endX)
		indices = np.where(mask)[0]
		if len(indices) == 0:
			raise ValueError("There are no data points in the range!")
		y_range=self.yValueList[mask]
		y_offset = np.mean(y_range)
		self.yValueList-=y_offset

	def eliminateInitialEffect(self,startX:float,endX:float):
		"""
		By eliminating the initial interval effects, such as friction influence, the data within the initial interval
		is essentially zeroed out.
        ----------------------------------------------------------------------------------------------------------------
        Inputs:
        		startX(float)-The start value of the range
        		endX(float)-The end value of the range
		"""
		mask = (self.xValueList >=startX) & (self.xValueList<=endX)
		if np.sum(mask) == 0:
			raise ValueError(f"There are no points in the range {startX}-{endX}!")
		friction_force = np.mean(np.abs(self.yValueList[mask]))
		###################
		dx = np.diff(self.xValueList, prepend=self.xValueList[0])
		dx = np.convolve(dx, np.ones(5) /5, mode='same')
		direction = np.sign(dx)
		# 填充为0的点（保持前一时刻方向）
		for i in range(1, len(direction)):
			if direction[i] == 0:
				direction[i] = direction[i - 1]
		if direction[0] == 0:
			direction[0] = 1
		#####################







########################################################################################################################
########################################################################################################################
if __name__=="__main__":
	testData = np.loadtxt('testSMACableData.txt')
	testDisp = list(testData[:, 0])
	testForce = list(testData[:, 1] * 1000)
	# plt.plot(testDisp, testForce)
	# plt.show()
	instance=HystereticCurveAnalysis(xDataList=testDisp,yDataList=testForce)
	instance.plotHystereticCurve(saveFig=False,multiColors=False)
	# instance.skeletonCurve(saveData=True,saveFig=False)
	# instance.plotLoop(loopNumber=2,saveData=False,saveFig=False)












