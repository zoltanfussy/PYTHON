#SQLite
#statements - snippets always end with ;
#clauses - commands are capitalized
SELECT * #all data
SELECT name AS 'Titles' #rename column, must use single quotes
SELECT DISTINCT tools #returns unique values from the column tools
FROM browse #from this database
LIMIT 10; #limit to 10 records
ORDER BY name DESC #by column "name" by descending order; default is ASC for ascending
#parameter in brackets ( ) defines input data, values can be passed instead of data type
CREATE TABLE table_name (
   id INTEGER, 
   column_2 data_type, 
   column_3 data_values
);
INSERT INTO table_name (id, column_2, column_3) #defines which columns will be altered
VALUES (1, 'Justin Bieber', 22); #one by one the values to a row
#modify table structure
ALTER TABLE table_name 
ADD COLUMN column_4 TEXT;
#change values
UPDATE table_name 
SET column_4 = '@taylorswift13' 
WHERE id = 4; #can use comparison, such as > < != >=

#alternatively:
WHERE name LIKE 'Se_en' # is able to return *similar* values with _ as wildcart for one character
WHERE name LIKE %man% #any or no characters may precede or follow, not case sensitive
WHERE name BETWEEN 'A' AND 'J' #not including values starting with J - because no name is just J
	AND value < 1000 #or OR
#remove rows
DELETE FROM table_name 
WHERE column_4 IS NULL;

#if-logic using CASE
SELECT name, depth,
	CASE
		WHEN size > 50 THEN 'Large stuff'
		WHEN size > 10 THEN 'Medium stuff'
		ELSE 'Small stuff'
	END AS 'Size container' #this effectively adds a column that evaluates the values covered
FROM data;

#constraints - insert only data of a given type
#SQL data types: INTEGER, TEXT, DATE, REAL
#PRIMARY KEY - can be used to uniquely identify rows, new rows with the same ID cannot be inserted
#UNIQUE - similar to PRIMARY KEY, only unique values will be retained; unlike PK, more columns can be UNIQUE
#NOT NULL - a value must be present
CREATE TABLE celebs (
   id INTEGER PRIMARY KEY, 
   name TEXT UNIQUE,
   date_of_birth TEXT NOT NULL,
   date_of_death TEXT DEFAULT 'Not Applicable'
);

COUNT()		# count the number of rows
SUM()		# the sum of the values in a column
MAX()/MIN() # the largest/smallest value
AVG()		# the average of the values in a column
ROUND(column, 0)		# round the values in the column to the number of decimals set
#eg. 
COUNT(*)
FROM table_name
WHERE column_3 = 0;

#calculate stuff for grouped data
SELECT category, COUNT(*) 
FROM fake_apps
GROUP BY category #or ROUND(price, 1)
ORDER BY category;
#you can also group by ROUND(rating) by using:
GROUP BY 1 				#1 refers to SELECTed! column number
HAVING COUNT(name) > 5	#much like WHERE, but filters groups instead of rows
						#HAVING comes after GROUP BY but before ORDER BY and LIMIT
#============
#an example code to calculate how many people enrolled and cancelled in a month - march
#CASE = if
SELECT COUNT(DISTINCT user_id) AS 'enrollments',
  COUNT(CASE
        WHEN strftime("%m", cancel_date) = '03'
        THEN user_id
  END) AS 'march_cancellations',
  ROUND(100.0 * COUNT(CASE
        WHEN strftime("%m", cancel_date) = '03'
        THEN user_id
  END) / COUNT(DISTINCT user_id)) AS 'churn_rate'
FROM pro_users
WHERE signup_date < '2017-04-01'
  AND (
    (cancel_date IS NULL) OR
    (cancel_date > '2017-03-01')
  );

#combining tables with identical number of columns:
SELECT *
FROM table1
UNION
SELECT *
FROM table2;

#combining tables by id:
SELECT *
FROM table1
JOIN table2 		#note this is a so-called inner join - only overlapping sets are joined
  ON table1.c2 = table2.c2;

LEFT JOIN table2	#this will only keep unshared values from table1
  ON table1.c2 = table2.c2
WHERE table2.c2 IS NULL;	#this will only keep unique values of table1


CROSS JOIN table2	#combinatorial join of two table columns, useful to subset items with a range of values
#an example code to calculate number of subscribers per each month based on start and end of their subscription:
SELECT month,
  COUNT(*)
FROM newspaper
CROSS JOIN months
WHERE start_month <= month AND end_month >= month
GROUP BY month;


WITH previous_query AS ( 	#calculate a temporary column to be used later
...
)

#example:
WITH previous_query AS (
	SELECT customer_id, COUNT(subscription_id) AS 'subscriptions'
	FROM orders
	GROUP BY customer_id
)
SELECT customers.customer_name, previous_query.subscriptions
FROM previous_query
JOIN customers
	ON previous_query.customer_id = customers.customer_id;