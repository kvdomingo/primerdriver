const path = require("path"),
    HtmlWebpackPlugin = require("html-webpack-plugin"),
    MiniCssExtractPlugin = require("mini-css-extract-plugin"),
    CleanWebpackPlugin = require("clean-webpack-plugin").CleanWebpackPlugin;


module.exports = {
    context: __dirname,
    entry: {
        main: ["./frontend/static/frontend/js/index"]
    },
    output: {
        path: path.resolve(__dirname, "frontend/static/frontend/bundles/"),
        publicPath: "/static/frontend/bundles/",
        filename: "main.[hash].js",
        chunkFilename: "[id].main.[hash].js",
        crossOriginLoading: "anonymous"
    },
    module: {
        rules: [
            {
                test: /\.(js|jsx)$/,
                exclude: [
                    /node_modules/,
                    /Content/,
                ],
                use: ["babel-loader"]
            },
            {
                test: /\.scss$/,
                use: ["style-loader", "css-loader", "sass-loader"]
            },
            {
                test: /\.css$/,
                use: [MiniCssExtractPlugin.loader, "css-loader"]
            },
            { test: /\.eot(\?v=\d+\.\d+\.\d+)?$/, loader: "file-loader" },
            { test: /\.(woff|woff2)$/, loader: "url-loader" },
            {
                test: /\.ttf(\?v=\d+\.\d+\.\d+)?$/,
                loader: "url-loader?mimetype=application/octet-stream"
            },
            {
                test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
                issuer: {
                    test: /\.jsx$/
                },
                use: ["babel-loader", "@svgr/webpack", "url-loader"]
            },
            {
                test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
                loader: "url-loader"
            },
            {
                test: /\.png(\?v=\d+\.\d+\.\d+)?$/,
                loader: "url-loader?limit=10000&mimetype=image/png"
            },
            {
                test: /\.gif(\?v=\d+\.\d+\.\d+)?$/,
                loader: "url-loader?limit=10000&mimetype=image/gif"
            }
        ]
    },
    plugins: [
        new CleanWebpackPlugin(),
        new MiniCssExtractPlugin({
            filename: "main.[hash].css",
            chunkFilename: "[id].main.[hash].css"
        }),
        new HtmlWebpackPlugin({
            template: path.resolve(__dirname, "frontend/jinja2/frontend/index.html.j2"),
            filename: "index.html"
        })
    ]
};
