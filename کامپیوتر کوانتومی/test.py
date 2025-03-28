import cirq

# 1. ایجاد کیوبیت‌ها
a = cirq.LineQubit(0)  # ورودی اول
b = cirq.LineQubit(1)  # ورودی دوم
carry = cirq.LineQubit(2)  # بیت نقلی (carry)
sum_bit = cirq.LineQubit(3)  # بیت جمع (sum)

# 2. ایجاد مدار کوانتومی
circuit = cirq.Circuit()

# 3. تنظیم مقادیر اولیه‌ی a و b
# اگر می‌خواهید a = 1 و b = 1 باشد، از گیت X استفاده کنید.
circuit.append([
    cirq.X(a),  # تنظیم a به 1
    cirq.X(b)   # تنظیم b به 1
])

# 4. اعمال گیت‌ها برای جمع باینری
circuit.append([
    cirq.CCX(a, b, carry),  # گیت Toffoli برای محاسبه carry
    cirq.CNOT(a, sum_bit),  # گیت CNOT برای محاسبه sum
    cirq.CNOT(b, sum_bit),
    cirq.measure(a, key='a'),  # اندازه‌گیری a
    cirq.measure(b, key='b'),  # اندازه‌گیری b
    cirq.measure(sum_bit, key='sum'),  # اندازه‌گیری بیت جمع
    cirq.measure(carry, key='carry')  # اندازه‌گیری بیت نقلی
])

# 5. چاپ مدار
print("مدار کوانتومی برای جمع دو عدد ۱ بیتی:")
print(circuit)

# 6. شبیه‌سازی مدار
simulator = cirq.Simulator()
result = simulator.run(circuit, repetitions=1000)

# 7. نمایش نتایج
print("\nنتایج اندازه‌گیری:")
print("a:", result.histogram(key='a'))
print("b:", result.histogram(key='b'))
print("Sum:", result.histogram(key='sum'))
print("Carry:", result.histogram(key='carry'))
