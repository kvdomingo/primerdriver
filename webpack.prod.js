const merge = require("webpack-merge"),
    common = require("./webpack.common.js"),
    CompressionPlugin = require("compression-webpack-plugin"),
    BundleTracker = require("webpack-bundle-tracker"),
    webpack = require("webpack"),
    BrotliPlugin = require('brotli-webpack-plugin');


module.exports = merge(common, {
    context: __dirname,
    mode: "production",
    plugins: [
        new BundleTracker({
            path: __dirname,
            filename: "webpack-stats.json"
        }),
        new webpack.optimize.SplitChunksPlugin({
            chunks: "async",
            cacheGroups: {
                vendors: {
                    test: /[\\/]node_modules[\\/]/,
                    priority: -10
                }
            }
        }),
        new webpack.optimize.ModuleConcatenationPlugin(),
        new CompressionPlugin({
            filename: "[path].gz[query]",
            algorithm: "gzip",
            test: /\.(js|css|html)$/,
            threshold: 10240,
            minRatio: 0.7,
        }),
        new BrotliPlugin({
            filename: "[path].br[query]",
            algorithm: "brotliCompress",
            test: /\.(js|css|html|svg)$/,
            threshold: 10240,
            minRatio: 0.7,
        })
    ]
});
