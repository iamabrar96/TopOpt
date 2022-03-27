'''
Validate plot of Iteration vs Compliance for Volume Fraction and Filter Radius
    1. Plot iteration vs compliance for a specific Volume Fraction and Filter Radius
    2. Plot iteration vs compliance but for 3 Volume Fraction
    3. Plot iteration vs compliance but for 3 different Filter Radii
'''
from parameters import Parameters
import matplotlib.pyplot as plt
from simp_optimizer import SimpOptimizer


# 1. Set the parameters to visualize the plot of iterations vs compliance
max_iterations = 200
volume_fraction = 0.3
filter_radius = 1.5

# 2. Set the volume fractions for plotting iterations vs compliance
vf_1 = 0.3
vf_2 = 0.4
vf_3 = 0.5
vf_4 = 0.7

# 3. Set the filter Radii for plotting iterations vs compliance
rmin_1 = 1.5
rmin_2 = 2
rmin_3 = 4
rmin_4 = 7



def plot_iterationvscompliance(max_iterations,volume_fraction,filter_radius):
    ''' Plot iteration vs compliance for a specific Volume Fraction and filter radius '''
    params = Parameters()
    params.volfrac = volume_fraction
    params.rmin = filter_radius
    params.max_loop = max_iterations
    fig,ax = plt.subplots()
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Volume Fraction = '+str(params.volfrac)+', Filter Radius = '+str(params.rmin))
    ax.set_title('Iteration vs Compliance')
    ax.set_xlabel('No of iterations')
    ax.set_ylabel('Compliance')
    ax.legend()
    plt.show()


def plot_iterationvscompliance_volfrac(max_iterations,vf_1,vf_2,vf_3,vf_4, filter_radius):
    '''Plot iteration vs compliance but for 3 different Volume Fractions'''
    params = Parameters()
    params.rmin = filter_radius
    params.max_loop = max_iterations
    params.volfrac = vf_1
    fig,ax = plt.subplots()
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Volume Fraction = '+str(params.volfrac))
    params.volfrac = vf_2
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Volume Fraction = '+str(params.volfrac))
    params.volfrac = vf_3
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Volume Fraction = '+str(params.volfrac))
    params.volfrac = vf_4
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Volume Fraction = '+str(params.volfrac))
    ax.set_title('Iteration vs Compliance')
    ax.set_xlabel('No of iterations')
    ax.set_ylabel('Compliance')
    ax.legend()
    plt.show()


def plot_iterationvscompliance_rmin(max_iterations,rmin_1,rmin_2,rmin_3, rmin_4, volume_fraction):
    '''Plot iteration vs compliance but for 3 different Filter Radii'''
    params = Parameters()
    params.rmin = rmin_1
    params.max_loop = max_iterations
    params.volfrac = volume_fraction
    fig,ax = plt.subplots()
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Filter Radius = '+str(params.rmin))
    params.rmin = rmin_2
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Filter Radius = '+str(params.rmin))
    params.rmin = rmin_3    
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Filter Radius = '+str(params.rmin))
    params.rmin = rmin_4
    op = SimpOptimizer(params)
    op.densityestimation()
    ax.plot(op.compliance,label = 'Filter Radius = '+str(params.rmin))
    ax.set_title('Iteration vs Compliance')
    ax.set_xlabel('No of iterations')
    ax.set_ylabel('Compliance')
    ax.legend()
    plt.show()


if __name__=='__main__':
    plot_iterationvscompliance(max_iterations,volume_fraction,filter_radius)
    plot_iterationvscompliance_volfrac(max_iterations,vf_1,vf_2,vf_3,vf_4,filter_radius)
    plot_iterationvscompliance_rmin(max_iterations,rmin_1,rmin_2,rmin_3, rmin_4, volume_fraction)