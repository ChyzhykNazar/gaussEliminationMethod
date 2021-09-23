# gaussEliminationMethod
📖Educational project, where presented solution for systems of linear equations using Gaussian elimination method

📖Навчальний проект, де представлено розв'язок СЛАР(системи лінійних алгебраїчних рівнянь) за допомогою метода Гауса

# Українська версія(Ukrainian version)

Метод Гауса є одним з найпоширеніших методів рішення СЛАР. У його основі
лежить ідея послідовного виключення невідомих, що приводить вихідну систему до
трикутного виду, у якому всі коефіцієнти нижче головної діагоналі дорівнюють нулю.
Існують різні обчислювальні схеми, що реалізують цей метод. Найбільше поширення
мають схеми з вибором головного елемента по рядку, по стовпцю, або по всій матриці.
Необхідно розв’язати СЛАР виду:

![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%231.png)

Класичний метод Гауса (метод виключення) грунтується на доведенні матриці
коефіцієнтів системи (7.1) до трикутного вигляду:

![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%232.png)

і складається з двох етапів: прямого ходу і зворотного підставлення. Етап прямого ходу
закінчується, коли одне з рівнянь системи стає рівнянням з одним невідомим. Далі,
здійснюючи зворотне підставлення, знаходять всі невідомі. З метою зменшення
обчислювальної похибки використовують метод Гауса з вибором головного елементу.
Під головним елементом будемо розуміти максимальний по модулю елемент матриці А,
обраний по заданій множині рядків та стовпців. Алгоритм методу полягає у наступному.

Спочатку знаходиться головний елемент матриці А і шляхом перестановки
рівнянь та стовпців він розміщується на місці елемента a11(перестановку стовпців
потрібно запам’ятовувати для вірного співставлення змінних з отриманими
значеннями). Потім ділимо перше рівняння на a11, після чого воно отримує вигляд:

![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%233.png)

Помножимо послідовно (7.2) на а21, а31, ..., аn1 та віднімемо відповідно з 2-го, ... ,
n-го рівняння системи. Отримаємо СЛАР A¹x = B¹:


![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%234.png)

Її елементи розраховані за формулами: 


![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%235.png)


![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%236.png)


![](https://github.com/ChyzhykNazar/gaussEliminationMethod/blob/8550ffcb6d9e933e26e5582a82e6dc60ec7e2f4c/images/%237.png)

# English version (Англійська версія)

https://en.wikipedia.org/wiki/Gaussian_elimination

https://www.dummies.com/education/math/calculus/how-to-use-gaussian-elimination-to-solve-systems-of-equations/
