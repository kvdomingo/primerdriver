const merge = require("webpack-merge"),
    common = require("./webpack.common.js"),
    BundleTracker = require("webpack-bundle-tracker");


module.exports = merge(common, {
    context: __dirname,
    mode: "development",
    devtool: "inline-source-map",
    plugins: [
        new BundleTracker({
            path: __dirname,
            filename: "webpack-stats.json",
            indent: 4
        })
    ]
});
