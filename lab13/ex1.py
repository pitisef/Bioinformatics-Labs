import tkinter as tk
import matplotlib.pyplot as plt

def multiply(matrix, vector):
    rows = len(matrix)
    cols = len(matrix[0])
    result = []
    
    for i in range(rows):
        total = 0
        for j in range(cols):
            total += matrix[i][j] * vector[j]
        result.append(total)
        
    return result

def run():
    initial_vector = [100, 50, 20]
    
    matrix = [
        [0.8, 0.1, 0.1],
        [0.2, 0.7, 0.1],
        [0.1, 0.3, 0.6]
    ]
    
    steps = 5
    current = initial_vector
    history = [current]
    
    for i in range(steps):
        next_val = multiply(matrix, current)
        history.append(next_val)
        current = next_val
        print(f"Step {i+1}: {current}")
        
    x_axis = range(len(history))
    
    plt.figure(figsize=(10, 6))
    
    comp1 = [v[0] for v in history]
    comp2 = [v[1] for v in history]
    comp3 = [v[2] for v in history]
    
    plt.plot(x_axis, comp1, label="Component 1")
    plt.plot(x_axis, comp2, label="Component 2")
    plt.plot(x_axis, comp3, label="Component 3")
    
    plt.xlabel("Step")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True)
    
root = tk.Tk()
root.title("Matrix Predictor")
root.geometry("800x600")

run()

root.mainloop()