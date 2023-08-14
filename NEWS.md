# construction 0.4.1

* Increase the default consensus threshold to 66%.

# construction 0.4.0

* The `construct()` function no longer uses the `threshold` argument. 
This argument is now used with the `consensus()` method that calculates consensus calculations based on a provided `adduct` and `threshold` arguments.

* The `construction()` method now supports parallel processing using the [`future`](https://future.futureverse.org/) package.

* The search logic within the `construction()` method has been streamlined to improve performance.

* The argument `library_path` in the method `construction()` now specifies the full path to the construction library.

* Added the `plotSankey()` method to plot structural overviews 

* Added the `structural_classifications` example data set for the `plotSankey()` plot example.

# construction 0.3.1

* Remove throttling of ClassyFire database queries.

* Added elapsed time to the progress bar.

# construction 0.3.0

* The package documentation is now available at <https://jasenfinch.github.io/construction/>.

* Function documentation improvements.

* Added support for the caching utilities now available in [`classyfireR`](https://aberhrml.github.io/classyfireR/) for caching structural classifications.

* Added a `construction()` method for the [`Assignments`](https://aberhrml.github.io/assignments/reference/Assignment-class.html) S4 class.

* Added the `Construction` S4 class, which is returned when using the `construction()` method for the [`Assignments`](https://aberhrml.github.io/assignments/reference/Assignment-class.html) S4 class.

* Added the molecular formula search number and percentage progress to the console output of `construction()`.

# construction 0.2.10

* Fixed error when 0 CIDs returned from Pubchem.

# construction 0.2.9

* Added a `NEWS.md` file to track changes to the package.

* Added Remotes field to DESCRIPTION.

* GitHub actions now used for continuous integration.
